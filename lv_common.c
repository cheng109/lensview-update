/**************************
LENSVIEW common routines for Grav Lensing
Copyright (C) Randall Wayth. 2006.
$Id: lv_common.c,v 1.34 2008/10/31 20:43:10 rwayth Exp rwayth $
***************************/
#include    <limits.h>
#include	<unistd.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>
#include	<signal.h>
#include    <float.h>
#include	"fitsio.h"
#include	"common.h"
#include	"log.h"
#include	"lv_geom.h"
#include	"lv_common.h"
#include	"lv_lens.h"
#include	"lv_image.h"

#define	POINT_MAP_CACHE_SIZE	4

/*************************
 * public global variables
 * ***********************/
int		g_ShutdownFlag =0;
float	g_FixedLambda =-1;
float	g_DataMean=0;
int     g_iDebugImgs=FALSE;
float   g_fTargetChiSqu=-1;

/***************************
Private function prototypes
***************************/
static	int	lv_mapPolygon(lv_lensmodel_t *pLens, lv_polygon *pImg1, lv_polygon *pSrc1, bool *pSkipVertex);
static  int	lv_updateMappingMatrix(lv_mapmatrix_t *pMap, lv_axissize_t ImgX, lv_axissize_t ImgY, lv_polygon *pSrcTri1);
static	int	lv_updateMappingWeight(lv_mapmatrix_t *pMap,lv_axissize_t imgX, lv_axissize_t imgY, lv_axissize_t srcX, lv_axissize_t srcY, real_t fWeight);
static	int	lv_searchCache(lv_point *pPoint);
static	void	lv_addToCache(lv_point *pImgPoint, lv_point *pSrcPoint, bool iRes);
static  void    lv_fpesighandle(int sig);
static	int		lv_CheckGridLensCentre(lv_lensmodel_t *pLens, lv_polygon *pImg1);

/***************************
Private global variables
****************************/
static	lv_point	imgPointCache[POINT_MAP_CACHE_SIZE];
static	lv_point	srcPointCache[POINT_MAP_CACHE_SIZE];
static	bool		srcSkipVertexCache[POINT_MAP_CACHE_SIZE];
static	short		iCacheLastEntry=-1,iCacheSize=0;

/***************************
****************************/

/***************************
Function:		lv_createMappingMatrix
Description:	Creates a 4-D mapping matrix between grid points in the source and image plane. Does not actually
				project any image. Only needs the axis and pixel angular sizes from pSourceImage and pProjImage, it doesn't
				matter if they have been formed yet.
Arguments:		pLens: lens model used to define mapping
				pMap:	resultant 4-D mapping matrix
				pSourceImage: source imgage details
				pProjImage: (yet to be) projected image details
Returns:		0- success. non-zero: failure
****************************/
int	lv_createMappingMatrix(lv_lensmodel_t *pLens, lv_mapmatrix_t *pMap, lv_image_t *pSourceImage, lv_image_t *pProjImage) {

	int		i,j,k,iResult=0, iStatus=0;
	real_t	x,y,x1,y1,fHalfImgRows,fHalfImgCols,fHalfSrcRows,fHalfSrcCols, fInvSrcPixelAngSize;
	lv_polygon	ImgGrid, SrcGrid, ImgTri1, ImgTri2, SrcTri1, SrcTri2, SrcPlane,tempTri1,tempTri2;
	bool	bSkipTri1=FALSE, bSkipTri2=FALSE, bSkipVertex[POINT_MAP_CACHE_SIZE+1];

	TRACE_IN(lv_createMappingMatrix);

	if (pSourceImage == NULL) {
		LOG_ERR("No source image. Exiting.");
		iStatus = 1;
		goto EXIT;
	}

	if (pProjImage == NULL) {
		LOG_ERR("No projected image. Exiting.");
		iStatus = 1;
		goto EXIT;
	}

	if (pProjImage->pMask == NULL) {
		pProjImage->pMask = (bool *) calloc(img_CalcImgSize(pProjImage),sizeof(bool));
		if (pProjImage->pMask==NULL) {
			LOG_ERR("Cannot create pixel mask for projected image");
			iStatus = 1;
			goto EXIT;
		}
	}

	fHalfImgRows = (real_t) (pProjImage->pAxisSize[1] / 2);
	fHalfImgCols = (real_t) (pProjImage->pAxisSize[0] / 2);
	fHalfSrcRows = (real_t) (pSourceImage->pAxisSize[1] / 2);
	fHalfSrcCols = (real_t) (pSourceImage->pAxisSize[0] / 2);

	/* since the pixel ang size doesn't change, this will save us doing lots
		of floating divisions in the loop. Do floating mult instead which is
		faster */
	if (pSourceImage->fPixelAngSize != 0) {
		fInvSrcPixelAngSize = 1/pSourceImage->fPixelAngSize;
	}
	else {
		TRACE(LOG_ERROR_PRI,"Source image pixel angular size is zero. Exiting...");
		goto EXIT;
	}

	/* make the polygons into triangles and rectangles */
	ImgTri1.numVertices = ImgTri2.numVertices = 3;
	SrcTri1.numVertices = SrcTri2.numVertices = 3;
	tempTri1.numVertices =tempTri2.numVertices = 3;
	SrcPlane.numVertices = ImgGrid.numVertices = SrcGrid.numVertices = 4;

	/* set up the src plane rectangle. 
	 * This is a big rectangle which is the size of the source plane
	 * used for checking if the projected source pixel lies entirely
	 * inside the source plane */
	SrcPlane.vertex[4].x = SrcPlane.vertex[1].x = (-pSourceImage->pAxisOffset[0]-fHalfSrcCols)*pSourceImage->fPixelAngSize;
	SrcPlane.vertex[2].y = SrcPlane.vertex[1].y = (-pSourceImage->pAxisOffset[1]-fHalfSrcRows)*pSourceImage->fPixelAngSize;
	SrcPlane.vertex[3].x = SrcPlane.vertex[2].x = (pSourceImage->pAxisSize[0]-fHalfSrcCols-pSourceImage->pAxisOffset[0])*pSourceImage->fPixelAngSize;
	SrcPlane.vertex[4].y = SrcPlane.vertex[3].y = (pSourceImage->pAxisSize[1]-fHalfSrcRows-pSourceImage->pAxisOffset[1])*pSourceImage->fPixelAngSize;
	geomSortPolygon(&SrcPlane);

	/* for each square in the image plane, working across then up */
	/* +ve Y is at top of image. +ve X is at right of image */
	for (j=0; j<pProjImage->pAxisSize[1]; j++) {

		sprintf(strMessage,"Starting row %d of %d",j,(int) (pProjImage->pAxisSize[1]-1));
		if (j == 0 || j == pProjImage->pAxisSize[1]-1) {
			TRACE(LOG_MED_PRI,strMessage);
		}
		else {
			TRACE(LOG_LOW_PRI,strMessage);
		}

		y = ( (real_t) j - fHalfImgRows ) * pProjImage->fPixelAngSize; 
		y1 = ( (real_t) j - fHalfImgRows + 1.0) * pProjImage->fPixelAngSize; 

		ImgTri1.vertex[2].y = ImgTri1.vertex[1].y = (real_t) j - fHalfImgRows;
		ImgTri1.vertex[3].y                       = ImgTri1.vertex[1].y + 1.0;

		ImgTri2.vertex[2].y = ImgTri2.vertex[3].y = ImgTri1.vertex[3].y;
		ImgTri2.vertex[1].y                       = ImgTri1.vertex[2].y;

		ImgGrid.vertex[1].y = y;
		ImgGrid.vertex[2].y = y;
		ImgGrid.vertex[3].y = y1;
		ImgGrid.vertex[4].y = y1;

		/* calculate the X and Y coords of the image pixel normalised to pixel units */
		for (i=0; i<pProjImage->pAxisSize[0] ; i++) {

			if(pProjImage->bExternalMask == FALSE) {
				pProjImage->pMask[i + j*pProjImage->pAxisSize[0]] = FALSE;
			}

			x = ((real_t) i - fHalfImgCols) * pProjImage->fPixelAngSize; 
			x1 = ((real_t) i - fHalfImgCols + 1.0) * pProjImage->fPixelAngSize; 

			ImgGrid.vertex[4].x = ImgGrid.vertex[1].x = x;
			ImgGrid.vertex[3].x = ImgGrid.vertex[2].x = x1;

			/* break the grid rectangle into two triangles. Keep the image plane trianges
				in grid units. The ImgGrid variable gets pixel units */
			ImgTri1.vertex[3].x = ImgTri1.vertex[1].x = (real_t) i - fHalfImgCols;
			ImgTri1.vertex[2].x                       = ImgTri1.vertex[1].x + 1.0;

			ImgTri2.vertex[2].x = ImgTri2.vertex[1].x = ImgTri1.vertex[2].x;
			ImgTri2.vertex[3].x                       = ImgTri1.vertex[1].x;

			for (k=0; k< POINT_MAP_CACHE_SIZE+1; k++) {
				bSkipVertex[k]=FALSE;
			}

			/* map the image grid into the source plane, still in normalised units */
			iResult = lv_mapPolygon(pLens, &ImgGrid, &SrcGrid, bSkipVertex);
			if (iResult != 0) {
			    iStatus = iResult;
				goto EXIT;
			}

			sprintf(strMessage,"Img in angle units.  [%g,%g],[%g,%g],[%g,%g],[%g,%g]"
							,ImgGrid.vertex[1].x, ImgGrid.vertex[1].y
							,ImgGrid.vertex[2].x, ImgGrid.vertex[2].y
							,ImgGrid.vertex[3].x, ImgGrid.vertex[3].y
							,ImgGrid.vertex[4].x, ImgGrid.vertex[4].y);
			TRACE(LOG_LOW_PRI,strMessage);
			sprintf(strMessage,"Srcs in angle units. [%g,%g],[%g,%g],[%g,%g],[%g,%g]"
							,SrcGrid.vertex[1].x, SrcGrid.vertex[1].y
							,SrcGrid.vertex[2].x, SrcGrid.vertex[2].y
							,SrcGrid.vertex[3].x, SrcGrid.vertex[3].y
							,SrcGrid.vertex[4].x, SrcGrid.vertex[4].y);
			TRACE(LOG_LOW_PRI,strMessage);

			sprintf(strMessage,"Img in grid units.  [%g,%g],[%g,%g],[%g,%g],[%g,%g]"
							,ImgTri1.vertex[1].x, ImgTri1.vertex[1].y
							,ImgTri1.vertex[2].x, ImgTri1.vertex[2].y
							,ImgTri2.vertex[2].x, ImgTri2.vertex[2].y
							,ImgTri2.vertex[3].x, ImgTri2.vertex[3].y);
			TRACE(LOG_LOW_PRI,strMessage);

			/* how to choose the direction of the diagonal in the img grid ?
			 * If the square is folded onto itself, then we can just define the
			 * diagonal in the other direction and it will no longer be folded.
			 * Having an img grid square unfolded/untwisted is the most accurate, so
			 * we should strive for this if possible */
			/* default diagonal goes from top left to bottom right of each pixel. Test
			 * if either the top right or bottom left corner of the grid maps into the
			 * trianlge which it is not part of. This indicates we have a fold */

			/* create temporary triangles in source plane which we can use to see
			 * if the triangle has been folded */
			tempTri1.vertex[1] = SrcGrid.vertex[1];
			tempTri1.vertex[2] = SrcGrid.vertex[2];
			tempTri1.vertex[3] = SrcGrid.vertex[4];
			tempTri2.vertex[1] = SrcGrid.vertex[2];
			tempTri2.vertex[2] = SrcGrid.vertex[3];
			tempTri2.vertex[3] = SrcGrid.vertex[4];
			if (geomInsidePolygon(&SrcGrid.vertex[3],&tempTri1) != 0 ||
				geomInsidePolygon(&SrcGrid.vertex[1],&tempTri2) != 0) {
				/* swap the diagonal in the image plane pixel which the triangles
				 * are defined on. This stops folding of triangles */
				SrcTri1.vertex[1].x = SrcTri2.vertex[1].x =	SrcGrid.vertex[1].x	* fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[0];
				SrcTri1.vertex[1].y = SrcTri2.vertex[1].y =	SrcGrid.vertex[1].y	* fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[1];
				SrcTri2.vertex[2].x = 						SrcGrid.vertex[2].x * fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[0];
				SrcTri2.vertex[2].y = 						SrcGrid.vertex[2].y * fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[1];
				SrcTri2.vertex[3].x = SrcTri1.vertex[2].x = SrcGrid.vertex[3].x * fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[0];
				SrcTri2.vertex[3].y = SrcTri1.vertex[2].y = SrcGrid.vertex[3].y * fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[1];
				SrcTri1.vertex[3].x =						SrcGrid.vertex[4].x * fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[0];
				SrcTri1.vertex[3].y =						SrcGrid.vertex[4].y * fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[1];
				sprintf(strMessage,"Srcs in grid units (alternate fold). [%g,%g],[%g,%g],[%g,%g],[%g,%g]"
						,SrcTri1.vertex[1].x, SrcTri1.vertex[1].y
						,SrcTri2.vertex[2].x, SrcTri2.vertex[2].y
						,SrcTri1.vertex[2].x, SrcTri1.vertex[2].y
						,SrcTri1.vertex[3].x, SrcTri1.vertex[3].y);
				TRACE(LOG_LOW_PRI,strMessage);
			}
			else {
				/* use "default" division of image plane cell where triangle runs
				 * from top right to bottom left */
				SrcTri1.vertex[1].x =						SrcGrid.vertex[1].x	* fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[0];
				SrcTri1.vertex[1].y =						SrcGrid.vertex[1].y	* fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[1];
				SrcTri2.vertex[1].x = SrcTri1.vertex[2].x = SrcGrid.vertex[2].x * fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[0];
				SrcTri2.vertex[1].y = SrcTri1.vertex[2].y = SrcGrid.vertex[2].y * fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[1];
				SrcTri2.vertex[3].x = SrcTri1.vertex[3].x = SrcGrid.vertex[4].x * fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[0];
				SrcTri2.vertex[3].y = SrcTri1.vertex[3].y = SrcGrid.vertex[4].y * fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[1];
				SrcTri2.vertex[2].x =						SrcGrid.vertex[3].x * fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[0];
				SrcTri2.vertex[2].y =						SrcGrid.vertex[3].y * fInvSrcPixelAngSize	+ pSourceImage->pAxisOffset[1];
				sprintf(strMessage,"Srcs in grid units (default fold). [%g,%g],[%g,%g],[%g,%g],[%g,%g]"
						,SrcTri1.vertex[1].x, SrcTri1.vertex[1].y
						,SrcTri1.vertex[2].x, SrcTri1.vertex[2].y
						,SrcTri2.vertex[2].x, SrcTri2.vertex[2].y
						,SrcTri2.vertex[3].x, SrcTri2.vertex[3].y);
				TRACE(LOG_LOW_PRI,strMessage);
			}

			/* add the half axis sizes back so that 0,0 is at the bottom edge of the grid */
			SrcTri1.vertex[1].x += fHalfSrcCols;
			SrcTri1.vertex[2].x += fHalfSrcCols;
			SrcTri1.vertex[3].x += fHalfSrcCols;
			SrcTri1.vertex[1].y += fHalfSrcRows;
			SrcTri1.vertex[2].y += fHalfSrcRows;
			SrcTri1.vertex[3].y += fHalfSrcRows;

			SrcTri2.vertex[1].x += fHalfSrcCols;
			SrcTri2.vertex[2].x += fHalfSrcCols;
			SrcTri2.vertex[3].x += fHalfSrcCols;
			SrcTri2.vertex[1].y += fHalfSrcRows;
			SrcTri2.vertex[2].y += fHalfSrcRows;
			SrcTri2.vertex[3].y += fHalfSrcRows;

			/* calculate if this image pixel is completely in the source plane.
			 * if not, we should set the mask flag for this pixel */
			bSkipTri1 = bSkipTri2 = FALSE;
			{
				lv_point	p1,p2,p3,p4;
				real_t	maxx,minx,maxy,miny;

				maxx = MAX(SrcGrid.vertex[1].x,SrcGrid.vertex[2].x);
				maxx = MAX(SrcGrid.vertex[3].x,maxx);
				maxx = MAX(SrcGrid.vertex[4].x,maxx);
				minx = MIN(SrcGrid.vertex[1].x,SrcGrid.vertex[2].x);
				minx = MIN(SrcGrid.vertex[3].x,minx);
				minx = MIN(SrcGrid.vertex[4].x,minx);
				maxy = MAX(SrcGrid.vertex[1].y,SrcGrid.vertex[2].y);
				maxy = MAX(SrcGrid.vertex[3].y,maxy);
				maxy = MAX(SrcGrid.vertex[4].y,maxy);
				miny = MIN(SrcGrid.vertex[1].y,SrcGrid.vertex[2].y);
				miny = MIN(SrcGrid.vertex[3].y,miny);
				miny = MIN(SrcGrid.vertex[4].y,miny);
				p1.x = minx;
				p1.y = miny;
				p2.x = maxx;
				p2.y = maxy;
				p3.x = maxx;
				p3.y = miny;
				p4.x = minx;
				p4.y = maxy;

				/* mask out any pixels which are not completely in the source plane */
				if ( (geomInsidePolygon(&p1,&SrcPlane)) &&
					(geomInsidePolygon(&p2,&SrcPlane)) &&
					(geomInsidePolygon(&p3,&SrcPlane)) &&
					(geomInsidePolygon(&p4,&SrcPlane)) ) {
						/* set the mask to TRUE, ie don't mask */
					if (pProjImage->bExternalMask == FALSE) pProjImage->pMask[i + j*pProjImage->pAxisSize[0]] = TRUE;
				}

				/* calculate weights of any pixels which have any part of them inside the source plane */
				/* it is necessary to make the distinction between masking (above) and weights. We mask
				 * pixels which are not completely inside the source plane otherwise they would be counted
					in the chi-squ which is not fair if the pixel cannot properly affect the outcome of
					the chi-squ. For the weights, it is necessary to calculate for any pixel which is
					even partially in the source plane because pixels which are on the edge of the source
					plane which have only a small fraction of one image pixel mapping to them, will have
					a very small total weight. When reverse-projecting, the value gets divided by this
					small weight and drives the edge pixels to stupidly large values. Hence we calculate
					the total weight for the source pixel even though much of the mapped image pixel
					is outside the source plane
				 */
				if ( !geomInsidePolygon(&p1,&SrcPlane) &&
					!geomInsidePolygon(&p2,&SrcPlane) && 
					!geomInsidePolygon(&p3,&SrcPlane) && 
					!geomInsidePolygon(&p4,&SrcPlane) ) {
					bSkipTri1 = bSkipTri2 = TRUE;
					sprintf(strMessage,"Img pixel (%d,%d) maps completely outside src plane. Skipping...",i,j);
					TRACE(LOG_LOW_PRI,strMessage);
				}
			}

			/* calculate which triangles should be skipped and calculate the area of the
			 * projected triangle, regardless of whether it maps to the source plane or
			 * not. We can then use it to construct the image plane magnification map
			 * later */
			pMap->pSumImg[i+j*pProjImage->pAxisSize[0]]=0;

			if (bSkipVertex[1]==TRUE || bSkipVertex[2]==TRUE || bSkipVertex[3]==TRUE|| bSkipVertex[4]==TRUE) {
				bSkipTri1 = bSkipTri2 = TRUE;
				sprintf(strMessage,"skip vertex 1 true (%d,%d,%d). skipping img: (%d,%d)\n",(int)bSkipVertex[1],(int)bSkipVertex[2],(int)bSkipVertex[4],i,j);
				TRACE(LOG_LOW_PRI,strMessage);
			}
			else {
				pMap->pSumImg[i+j*pProjImage->pAxisSize[0]] += geomPolygonArea(&SrcTri1);
				pMap->pSumImg[i+j*pProjImage->pAxisSize[0]] += geomPolygonArea(&SrcTri2);
				sprintf(strMessage,"Img pixel (%d,%d) has total area in src plane %g",i,j,pMap->pSumImg[i+j*pProjImage->pAxisSize[0]]);
				TRACE(LOG_LOW_PRI,strMessage);
			}

			/* calculate the area of overlap */
			if (bSkipTri1 == FALSE) {
				iResult = lv_updateMappingMatrix(pMap,i,j, &SrcTri1);
			}
			if (bSkipTri2 == FALSE) {
				iResult = lv_updateMappingMatrix(pMap,i,j, &SrcTri2);
			}

			if (iResult != 0) {
				sprintf(strMessage,"lv_updateMappingMatrix failed for image co-ords (x,y) = (%d,%d)",i,j);
				TRACE(LOG_HIGH_PRI,strMessage);
				iStatus = iResult;
				goto EXIT;
			}
		}
	}
	/* calculate the % of the I,J weights which mapped to K,L values */
	/* calculate also the sums of areas for each img and source pixel 
	 * these will be used often in later calculations */

	TRACE(LOG_LOW_PRI,"Calculating fraction of mapping matrix used and array totals");
	
	for (i=0,j=0; i<pMap->dimension[0]*pMap->dimension[1]; i++) {
		int	k;

		if (pMap->array[i] != NULL) {
			j++;
			for(k=pMap->min[i];k <= pMap->max[i]; k++) {
                                pMap->pSumSrc[k] += pMap->array[i][k];
			}
		}
	}

	sprintf(strMessage,"%d percent of %d element I,J matrix was used. Total K,L memory used: %d K",j*100/(i+1),i,(int)(j*pMap->dimension[2]*pMap->dimension[3]*sizeof(real_t)/1024));
	TRACE(LOG_HIGH_PRI,strMessage);

EXIT:
	TRACE_OUT;
	return iStatus;
}

/***************************
Function:	lv_updateMappingMatrix
Description:	Calculate the area of overlap (in the source plane) 
				for a triangle projected from the image plane.
Arguments:	pMap: map/weights matrix.
			ImgX,ImgY: (x,y) co-ord of the image plane grid square
			pSrcTri: projected triangle
Returns:	0: success, non-zero: failure
****************************/
static	int	lv_updateMappingMatrix(lv_mapmatrix_t *pMap, lv_axissize_t ImgX, lv_axissize_t ImgY, lv_polygon *pSrcTri) {

	int				iResult=0;
	real_t			fXMin, fXMax, fYMin, fYMax, fArea=0;
	lv_axissize_t	iGridX, iGridY;
	lv_polygon		grid;

	TRACE_IN(lv_updateMappingMatrix);

	grid.numVertices = 4;

	/* find the X and Y bounds of the triangle in the source plane */
	fXMin = MIN(pSrcTri->vertex[1].x,pSrcTri->vertex[2].x);
	fXMin = MIN(fXMin,pSrcTri->vertex[3].x);
	fXMax = MAX(pSrcTri->vertex[1].x,pSrcTri->vertex[2].x);
	fXMax = MAX(fXMax,pSrcTri->vertex[3].x);

	fYMin = MIN(pSrcTri->vertex[1].y,pSrcTri->vertex[2].y);
	fYMin = MIN(fYMin,pSrcTri->vertex[3].y);
	fYMax = MAX(pSrcTri->vertex[1].y,pSrcTri->vertex[2].y);
	fYMax = MAX(fYMax,pSrcTri->vertex[3].y);

	/* make sure the bounds are actually in the source grid somewhere */
	fXMin = MAX(fXMin,0);
	fYMin = MAX(fYMin,0);
	fXMax = MIN(fXMax,(pMap->dimension[2]-1));
	fYMax = MIN(fYMax,(pMap->dimension[3]-1));

	sprintf(strMessage,"X min/max bounds: %g, %g. Y min/max bounds: %g, %g.",fXMin, fXMax, fYMin, fYMax);
	TRACE(LOG_LOW_PRI,strMessage);

	sprintf(strMessage,"Looking for grid overlap for img grid (%ld,%ld). Source plane triangle (%f,%f),(%f,%f),(%f,%f).", 
		ImgX,ImgY,
		pSrcTri->vertex[1].x, pSrcTri->vertex[1].y, pSrcTri->vertex[2].x, pSrcTri->vertex[2].y,
		pSrcTri->vertex[3].x, pSrcTri->vertex[3].y);
	TRACE(LOG_LOW_PRI,strMessage);

	/* for each grid square which COULD be overlapped by this,
		find the area of overlap */
	for (iGridX = (lv_axissize_t) floor(fXMin); iGridX <= (lv_axissize_t) floor(fXMax); iGridX++) {
		grid.vertex[1].x = iGridX;
		grid.vertex[3].x = grid.vertex[2].x = iGridX+1;
		grid.vertex[4].x = iGridX;

		for (iGridY = (lv_axissize_t) floor(fYMin); iGridY <= (lv_axissize_t) floor(fYMax); iGridY++) {
			fArea =0;
			grid.vertex[1].y = iGridY;
			grid.vertex[2].y = iGridY;
			grid.vertex[4].y = grid.vertex[3].y = iGridY+1;

			fArea = geomAreaOverlap(&grid,pSrcTri);

			sprintf(strMessage,"Calculating overlap for img grid (%ld,%ld), source grid square (%ld,%ld). Area: %g",ImgX,ImgY,iGridX,iGridY,fArea);
			TRACE(LOG_LOW_PRI,strMessage);

			if (fArea != 0) {
				/* update the mapping matrix */
				iResult = lv_updateMappingWeight(pMap, ImgX, ImgY, iGridX, iGridY, fArea);
			}
			if (iResult != 0) {
				goto EXIT;
			}
		}
	}
EXIT:
	TRACE_OUT;
	return iResult;
}

/***************************
Function:	lv_updateMappingWeight
Description:	Update a weight for a pixel in the source plane which has
				some part of an image plane trianlge mapped to it.
Arguments:	pMap: map/weights matrix,
			imgX,imgY,srcX,srcY: (x,y) co-ords of pixel in image and source planes
			fWeight: overlap area.
Returns:	0: success, non-zero: failure
****************************/
static	int	lv_updateMappingWeight(lv_mapmatrix_t *pMap, lv_axissize_t imgX, lv_axissize_t imgY, lv_axissize_t srcX, lv_axissize_t srcY, real_t fWeight) {
	mapwt_t	**pSrcArray;
	int		iResult =0;
	lv_axissize_t	iSrcOffset=0, iImgOffset =0,iNumElements=0;

	iSrcOffset = (srcY*pMap->dimension[2] + srcX);
	iImgOffset = (imgY*pMap->dimension[0] + imgX);
	pSrcArray = pMap->array + iImgOffset;
	if (*pSrcArray == NULL) {
/*
		fprintf(stderr,"Allocating source grid weights array for image grid %d, %d\n",(int) imgX,(int) imgY);
*/
		/* no weights have been assigned here so far */
		iNumElements = pMap->dimension[2] * pMap->dimension[3];
		if (iNumElements <= 0) {
			fprintf(stderr,"lv_updateMappingWeight Error: Attempt to allocate array of size (%ld)\n",(long) iNumElements);
			iResult = 1;
			goto EXIT;
		}

		*pSrcArray = (mapwt_t *) calloc(iNumElements , sizeof(mapwt_t));
		if (*pSrcArray == NULL) {
			fprintf(stderr,"lv_updateMappingWeight: malloc failed for %dx%d array of real_t for img grid %d, %d\n",(int) pMap->dimension[2],(int) pMap->dimension[3], (int) imgX, (int) imgY);
			iResult = 1;
			goto EXIT;
		}
	}
	*(*pSrcArray + iSrcOffset) += fWeight;
	if (iSrcOffset < pMap->min[iImgOffset]) {
		pMap->min[iImgOffset] = iSrcOffset;
	}
	if (iSrcOffset > pMap->max[iImgOffset]) {
		pMap->max[iImgOffset] = iSrcOffset;
	}
EXIT:
	return iResult;
}

/***************************
Function:	lv_mapPolygon
Description:	project a polygon's vertices to the source plane. Recalls the last few vertices
				which have been projected since each triangle in the image plane will have
				vertices in common with other triangles. Should reduce the overall number of
				"projections" by about 1/2.
Arguments:	pLens: lens model structure
			pImg1: polygon in image plane to project
			pSrc1: projected polygon in source plane
			pSkipVertex: flag to skip the vertex if it is at a special point
					like the centre of a singular lens
Returns:	0: success, non-zero: failure
****************************/
static	int	lv_mapPolygon(lv_lensmodel_t *pLens, lv_polygon *pImg1, lv_polygon *pSrc1, bool *pSkipVertex) {

	int		i,iCacheIndex, iStatus=0, bSkipGrid=FALSE;
	real_t	fMag;

/*
*	TRACE_IN(lv_mapPolygon);
*/

	/* check to see if this img plane grid contains a lens centre, 
	 * if so it needs special treatment */
	if (lv_CheckGridLensCentre(pLens,pImg1) != 0) bSkipGrid=TRUE;

	/* project into the source plane */
	for (i=1; i<=pImg1->numVertices; i++) {
		iCacheIndex = lv_searchCache(&(pImg1->vertex[i]));
		if (iCacheIndex < 0) {
			iStatus = lv_projectRay(pLens, pImg1->vertex[i].x, pImg1->vertex[i].y , &(pSrc1->vertex[i].x), &(pSrc1->vertex[i].y), &fMag);
            if (iStatus == LM_BADPROJ) return iStatus;
			lv_addToCache(&(pImg1->vertex[i]),&(pSrc1->vertex[i]),(iStatus == 0? FALSE: TRUE));
			pSkipVertex[i] = (iStatus==LM_IGNORE_PROJ);
			if (iStatus == LM_IGNORE_PROJ){
				iStatus=0;
			}
/*
*			TRACE(LOG_LOW_PRI,"Cache miss");
*/
		}
		else {
			pSrc1->vertex[i] = srcPointCache[iCacheIndex];
			pSkipVertex[i] = srcSkipVertexCache[iCacheIndex];
/*
*			TRACE(LOG_LOW_PRI,"Cache hit");
*/
		}
	}
	if (bSkipGrid) {
		for (i=0;i<=pImg1->numVertices; i++) pSkipVertex[i]=TRUE;
	}
/*
*	TRACE_OUT;
*/
	return iStatus;
}


/***************************
Function:	lv_CheckGridLensCentre
Description:	check to see if the pixel contains the centre of
				any singular lens component
Arguments:	pLens: lens model, pImg1: polygon defining the region of the pixel
Returns:	1- contains a lens centre, 0- doesn't
****************************/
static int	lv_CheckGridLensCentre(lv_lensmodel_t *pLens, lv_polygon *pImg1) {
	int	iStatus=0,i;
	lv_point p1;

	for (i=0; i< pLens->iNumComponents; i++) {
		p1.x = -g_PixelResn*(pLens->pComponent[i].fXoffset);
		p1.y = -g_PixelResn*(pLens->pComponent[i].fYoffset);
		switch(pLens->pComponent[i].iType) {
			case LM_SPEMD:
					if (pLens->pComponent[i].fParameter[4] > 0) break;
			case LM_PIEP:
					if (pLens->pComponent[i].fParameter[3] > 0) break;
			case LM_SIE:
			case LM_PTMASS:
			case LM_ISO_SPHERE:
					iStatus = geomInsidePolygon(&p1,pImg1);
					break;
			default:
				/* other lens models we don't care about */
				break;
		}
	}

	return iStatus;
}

/***************************
Function:	lv_searchCache
Description:	Search the cache of recently projected points and return
				the index if successful
Arguments:	pPoint: a point in the image plane to project
Returns:	index in cache or -1 on failure
****************************/
static	int	lv_searchCache(lv_point *pPoint) {
	register int i;
	int	iResult=-1;

	for(i=0; i<iCacheSize; i++) {
		if(pPoint->x == imgPointCache[i].x && pPoint->y == imgPointCache[i].y) {
			iResult = i;
			break;
		}
	}
	return iResult;
}

/***************************
Function:	lv_addToCache
Description:	Add a point to the projection cache
Arguments:	pImgPoint,pSrcPoint: point in image plane and source plane
			iRes: flag to skip point
Returns:	void
****************************/
static	void	lv_addToCache(lv_point *pImgPoint, lv_point *pSrcPoint, bool iRes) {

	iCacheLastEntry = (iCacheLastEntry+1) % POINT_MAP_CACHE_SIZE;
	imgPointCache[iCacheLastEntry] = *pImgPoint;
	srcPointCache[iCacheLastEntry] = *pSrcPoint;
	srcSkipVertexCache[iCacheLastEntry]= iRes;
	if (++iCacheSize > POINT_MAP_CACHE_SIZE) {
		iCacheSize = POINT_MAP_CACHE_SIZE;
	}
}


/***************************
Function:	lv_allocMappingMatrix
Description:	make an image mapping weights matrix. This is an abstraction in case I need to
				change implementation later
Arguments:	i,j,k,l: The axis sizes of the image (i,j) and source(k,l) arrays respectively
Returns:	pointer or NULL on failure
****************************/
lv_mapmatrix_t	*lv_allocMappingMatrix(lv_axissize_t ImgX,lv_axissize_t ImgY,lv_axissize_t SrcX,lv_axissize_t SrcY) {
	lv_mapmatrix_t	*pTemp=NULL;
	lv_axissize_t p;

	TRACE_IN(lv_allocMappingMatrix);

	sprintf(strMessage,"Allocating mapping matrix array for matrix sizes i: %d, j: %d, k: %d, l: %d",(int)ImgX,(int)ImgY,(int)SrcX,(int)SrcY);
	TRACE(LOG_MED_PRI,strMessage);

	if (ImgX*ImgY == 0 || SrcX*SrcY == 0) {
		LOG_ERR("Array dimension is zero. Exiting.");
		goto EXIT;
	}

	pTemp = (lv_mapmatrix_t *) malloc(sizeof(lv_mapmatrix_t));
	if (pTemp == NULL) {
		LOG_ERR("Failed to malloc lv_mapmatrix_t struct");
		goto EXIT;
	}

	/* change this ordering at your peril... */
	pTemp->dimension[0] = ImgX;
	pTemp->dimension[1] = ImgY;
	pTemp->dimension[2] = SrcX;
	pTemp->dimension[3] = SrcY;

	pTemp->array = (mapwt_t **) malloc (ImgX * ImgY*sizeof(mapwt_t *));
	if (pTemp->array == NULL) {
		sprintf(strMessage,"Failed to malloc weights array. Size %ld bytes",(long) sizeof(real_t *) * ImgX * ImgY );
		LOG_ERR(strMessage);
		free(pTemp);
		pTemp=NULL;
		goto EXIT;
	}
	sprintf(strMessage,"Allocated weights i,j array for %ld elements",(long) ImgX * ImgY);
	TRACE(LOG_MED_PRI,strMessage);

	pTemp->min = calloc(ImgX * ImgY,sizeof(lv_axissize_t));
	pTemp->max = calloc(ImgX * ImgY,sizeof(lv_axissize_t));
	pTemp->pSumImg = calloc(ImgX * ImgY,sizeof(real_t));
	pTemp->pSumSrc = calloc(SrcX * SrcY,sizeof(real_t));
	pTemp->pMultImgPix = calloc(ImgX * ImgY,sizeof(bool));

	if (pTemp->min == NULL || pTemp->max == NULL || pTemp->pSumImg==NULL
			|| pTemp->pSumSrc==NULL || pTemp->pMultImgPix==NULL) {
		LOG_ERR("Failed to calloc extra image space.");
		free(pTemp->array);
		free(pTemp);
		pTemp=NULL;
		goto EXIT;
	}

	/* init everything */
	for (p=0; p < ImgX * ImgY; p++) {
		pTemp->array[p] = NULL;
		pTemp->min[p] = INT_MAX;
		pTemp->max[p] = -1;
	}

EXIT:
	TRACE_OUT;
	return pTemp;
}

/***************************
Function:	lv_freeMappingMatrix
Description:	Free the memory for the weights matrix
Arguments:	pMap: pointer to the weights
Returns:	0: success
****************************/
int	lv_freeMappingMatrix(lv_mapmatrix_t *pMap) {

	int i, count=0;

	TRACE_IN(lv_freeMappingMatrix);

	if (pMap == NULL || pMap->array == NULL) {
		goto EXIT;
	}

	for (i=0; i< pMap->dimension[0] * pMap->dimension[1]; i++) {
		if (pMap->array[i] != NULL) {
			free(pMap->array[i]);
			count++;
		}
	}

	sprintf(strMessage,"Freed %d k,l arrays of size %d elements",count,(int) (pMap->dimension[2]*pMap->dimension[3]));
	TRACE(LOG_MED_PRI,strMessage);

	if (pMap->array != NULL) {
		free(pMap->array);
	}
	if (pMap->min != NULL) {
		free(pMap->min);
	}
	if (pMap->max != NULL) {
		free(pMap->max);
	}
	if (pMap->pSumImg != NULL) {
		free(pMap->pSumImg);
	}
	if (pMap->pSumSrc != NULL) {
		free(pMap->pSumSrc);
	}
	if (pMap->pMultImgPix !=NULL) {
		free(pMap->pMultImgPix);
	}
	free(pMap);

EXIT:

	TRACE_OUT;
	return 0;
}


/***************************
Function:	lv_zeroMapMatrix
Description: zero all the weights in the mapping matrix
Arguments:
Returns:
****************************/
int	lv_zeroMapMatrix(lv_mapmatrix_t *pMap) {
	lv_axissize_t			i,j,iSizeImg,iSizeSrcArr;

	TRACE_IN(lv_zeroMapMatrix);

	if (pMap == NULL) {
		goto EXIT;
	}

	TRACE(LOG_MED_PRI,"Zeroing map matrix.");

	iSizeImg = pMap->dimension[0] * pMap->dimension[1];
	iSizeSrcArr = pMap->dimension[2] * pMap->dimension[3];

	for(i=0; i<iSizeImg; i++) {
		if (pMap->array[i] != NULL) {
			for (j=0; j< iSizeSrcArr; j++) pMap->array[i][j] =0;
		}
		pMap->min[i] = INT_MAX;
		pMap->max[i] = -1;
	}
	lv_ZeroRealArray(pMap->pSumImg,iSizeImg);
	lv_ZeroRealArray(pMap->pSumSrc,iSizeSrcArr);
	for(i=0; i<iSizeImg; i++) pMap->pMultImgPix[i]=FALSE;

EXIT:
	TRACE_OUT;
	return 0;
}


/***************************
Lens Ray Tracing
****************************/

/***************************
Function:	lv_projectSourceThruMapMatrix
Description: projects a given source image into the image struct using the
			(already calculated) mapping weights.
Arguments:	pSource- source image to be projected
			pImage-	destination for projection. Must already be allocated etc.
			pMap-	mapping matrix weights.
Returns:	0: success, nonzero: failure.
****************************/
int	lv_projectSourceThruMapMatrix(lv_image_t *pSource, lv_image_t *pImg, lv_mapmatrix_t *pMap) {

	register	lv_axissize_t	j;
	register	real_t	fWeightedSumKL, fWeight;
	mapwt_t		*pSourceWeight;
	lv_axissize_t		i, iSizeImgArr = 0, iSizeSrcArr = 0;
	real_t		*pImgArray,*pSourceVal;
	int			iStatus =0;

	TRACE_IN(lv_projectSourceThruMapMatrix);

	/* for each i,j grid, take the weighted sum of source vals 
		divided by the sum of the weights */

	pImgArray = pImg->pImage;
	iSizeImgArr = pMap->dimension[1] * pMap->dimension[0];
	iSizeSrcArr = pMap->dimension[2] * pMap->dimension[3];

	if (iSizeImgArr == 0 || iSizeSrcArr == 0) {
		sprintf(strMessage,"Size of source/img array is zero. Exiting. Source: %d, image: %d", (int) iSizeSrcArr, (int) iSizeImgArr);
		TRACE(LOG_HIGH_PRI,strMessage);
		iStatus = 1;
	}

	/* loop through image array */
	for (i=0; i< iSizeImgArr; i++) {	

		/* if this pixel cannot be affected by anything in the source, then we should mark it so */
		pImgArray[i] = NULL_PIXEL;
		fWeightedSumKL =0.0;

		/* project the weighted source for this image pixel position */
		pSourceWeight = pMap->array[i];

		if (pSourceWeight != NULL) {

			pSourceVal = pSource->pImage;

			for(j=pMap->min[i]; j<=pMap->max[i]; j++) {
				fWeight = pSourceWeight[j];
				if (fWeight != 0.0) {
					fWeightedSumKL += fWeight * (pSourceVal[j] - pSource->fDefaultVal);
				}
			}

		}

		if (pMap->pSumImg[i] != 0.0) {
			pImgArray[i] = (real_t) fWeightedSumKL / pMap->pSumImg[i];
		}
	}


	TRACE_OUT;
	return iStatus;
}

/***************************
Function:	lv_projectRay
Description:	Project a single ray from source to image
Arguments:		pLens: Lens model
				pImg: Image plane data. Needed for angular size of pixels in image plane
				fImgX, fImgY: x and y coords of the pixel in the image plane in angular units
				pSrcY, pSrcX: result projected location in the source plane
				pMag: magnification at that point. (not currently implemented)
Returns:	0:	success, non-zero: failure
****************************/
int	lv_projectRay(lv_lensmodel_t *pLens, real_t fImgX, real_t fImgY, real_t *pSrcX, real_t *pSrcY, real_t *pMag){

	int	iStatus =0;
	real_t	fDeflectionX=0.0, fDeflectionY = 0.0;

/*
	TRACE_IN(lv_projectRay);
*/

	*pSrcX = 0; *pSrcY = 0;

	/* calculate the deflection angle vector */
	/* the defleciton is already the scaled deflection angle. */
	iStatus = lm_CalcDeflection(pLens, fImgX, fImgY, &fDeflectionX, &fDeflectionY, pMag);
	if (iStatus != 0 && iStatus != LM_IGNORE_PROJ) goto EXIT;

	/* add deflection to source position */
	/* remember that deflection is in -R direction, so make direction of X/Y negative */
	*pSrcX = fImgX - fDeflectionX;
	*pSrcY = fImgY - fDeflectionY;

/*
	sprintf(strMessage,"Img x,y: %g, %g. Source x,y: %g, %g",fImgX, fImgY, *pSrcX, *pSrcY);
	TRACE(LOG_LOW_PRI,strMessage);
	TRACE_OUT;
*/
EXIT:
	return iStatus;
}

/***************************
Function:	lv_projectImage
Description:	project a source image onto the image plane using the lens model
				provided using straight ray-tracing
Arguments:	pLens: pointer to lens model
			pSourceImage: pointer to source image which will be projected
			pProjImage:	pointer to projected image array to be created
Returns:	0:	success
****************************/
int	lv_projectImage(lv_lensmodel_t *pLens, lv_image_t *pSourceImage, lv_image_t *pProjImage) {

	int		iStatus =0,i,j,iNumRows =0, iNumCols =0, SourceI, SourceJ;
	real_t	fNormX=0, fNormY=0, fMagnification = 1, fSourceX, fSourceY;
	real_t	*pSource=NULL, *pProj=NULL;

	TRACE_IN(projectImage);

	if (pSourceImage == NULL) {
		TRACE(LOG_HIGH_PRI,"argument pointer pSourceImage is NULL. Exiting.")
		iStatus =1;
		goto EXIT;
	}

	if (pProjImage == NULL) {
		TRACE(LOG_HIGH_PRI,"argument pointer pProjImage is NULL. Exiting.")
		iStatus =1;
		goto EXIT;
	}

	iNumRows = pProjImage->pAxisSize[1];
	iNumCols = pProjImage->pAxisSize[0];

	/* check the pixel angular size vs the angular size of the Einstein Radius (ER) */
	sprintf(strMessage,"Ang size of image pixel is %g, total ang sizeof image is %g, ang size of source is %g", (double) pProjImage->fPixelAngSize, (double) (pProjImage->fPixelAngSize * iNumCols), (double) pSourceImage->fPixelAngSize);
	TRACE(LOG_MED_PRI,strMessage);

	/* work along rows then through the rows in the image. */

	for (j=0; j < iNumRows; j++) {
		fNormY =  (j - iNumRows/2) * pProjImage->fPixelAngSize; 
		/* work across the row */
		for (i=0; i < iNumCols; i++) {

			/* calculate the X and Y coords of the image pixel normalised to pixel angle units */
			fNormX =  (i - iNumCols/2) * pProjImage->fPixelAngSize; 

			/* project the Ray from Image to source using angle units*/
			lv_projectRay(pLens, fNormX, fNormY, &fSourceX, &fSourceY, &fMagnification);

			/* calculate the source grid positions. This is with 0,0 in the centre of the grid */
			SourceI = fSourceX/pSourceImage->fPixelAngSize;
			SourceJ = fSourceY/pSourceImage->fPixelAngSize;

			sprintf(strMessage,"Angle units. Img x,y: %g, %g.\tSrc x,y: %g,%g",fNormX,fNormY,fSourceX,fSourceY);
			TRACE(LOG_LOW_PRI,strMessage);

			sprintf(strMessage,"Grid units.  Img x,y: %d, %d.\tSrc x,y: %d,%d",i-iNumCols/2,j-iNumRows/2,SourceI,SourceJ);
			TRACE(LOG_LOW_PRI,strMessage);

			/* add the offsets of the grid length to put 0,0 at bottom left */
			SourceI += pSourceImage->pAxisSize[0]/2;
			SourceJ += pSourceImage->pAxisSize[1]/2;

			/* if the X,Y coords of the source pixel are actaully in the image, copy the value into the projected image */
			if (SourceI >= 0 && SourceI <pSourceImage->pAxisSize[0] && SourceJ >=0 && SourceJ < pSourceImage->pAxisSize[1] ) {

				/* point to the start of the image */
				pSource = pSourceImage->pImage;
				pProj = pProjImage->pImage;

				/* calculate the offset */
				pSource += (SourceJ* pSourceImage->pAxisSize[0]) + SourceI;
				pProj += (iNumRows * j) + i;
				*pProj = *pSource * fMagnification;
			}
		}
	}

EXIT:
	TRACE_OUT;
	return iStatus;
}

/**********************
File management
***********************/

/***************************
Function:	lv_write_image_to_file
Description:	writes a lv_image_t struct to a FITS file.
Arguments:	pImagePtr:	image struct to be written
			strFileName: file name to be written to.
			bClobber: overwrite existing file if it exists
Returns:	0: success, nonzero: failure
****************************/
int	lv_write_image_to_file(lv_image_t *pImagePtr, char	*strFileName, bool	bClobber) {
	fitsfile	*pFitsFile = NULL;
	FILE		*fp=NULL;
	int	iStatus=0;

	TRACE_IN(lv_write_image_to_file);
	if (strFileName == NULL || pImagePtr == NULL) {
		sprintf(strMessage, "Input pointer(s) are NULL. Image: %p, Filename: %p",pImagePtr,strFileName);
		TRACE(LOG_HIGH_PRI,strMessage);
		iStatus =1;
		goto EXIT;
	}
	if ((fp = fopen(strFileName, "r")) != NULL) {
		fclose(fp);

		if (bClobber == TRUE) {
			sprintf(strMessage,"Overwriting image file <%s>",strFileName);
			TRACE(LOG_MED_PRI,strMessage);
			remove(strFileName);
		}
		else {
			sprintf(strMessage,"File <%s> already exists. Leaving it.",strFileName);
			TRACE(LOG_MED_PRI,strMessage);
			goto EXIT;
		}
	}
	else {
		sprintf(strMessage,"Writing image to file name <%s>",strFileName);
		TRACE(LOG_MED_PRI,strMessage);
	}

	fits_create_file(&pFitsFile, strFileName, &iStatus);
	if (iStatus != 0) {
		sprintf(strMessage,"fits_create_file failed with status %d",iStatus);
		TRACE(LOG_HIGH_PRI,strMessage);
		goto EXIT;
	}

	fits_create_img(pFitsFile, pImagePtr->iFitsBitPix, pImagePtr->iNumAxes, pImagePtr->pAxisSize, &iStatus);
	if (iStatus != 0) {
		sprintf(strMessage,"fits_create_img failed with status %d",iStatus);
		TRACE(LOG_HIGH_PRI,strMessage);
		goto FILE_OPEN_EXIT;
	}

	fits_write_img(pFitsFile, pImagePtr->iTableDataType, 1, pImagePtr->pAxisSize[0]*pImagePtr->pAxisSize[1], pImagePtr->pImage, &iStatus);
	if (iStatus != 0) {
		sprintf(strMessage,"fits_write_img failed with status %d",iStatus);
		TRACE(LOG_HIGH_PRI,strMessage);
	}

FILE_OPEN_EXIT:
	fits_close_file(pFitsFile,&iStatus);

EXIT:
	TRACE_OUT;
	return iStatus;
}

/***************************
 Added by Jun Cheng
 Function:	lv_residual_img_file
 Description:	Get residual of two fits images
 Arguments:	 <lv_image_t input> <lv_image_t model>
 Returns:	<lv_image_t residual> or NULL on failure.
 ****************************/
lv_image_t * lv_residual_img_file(lv_image_t *inputImage, lv_image_t *modelImage) {
    lv_image_t * residulImage = lv_duplicate_image(inputImage);
/*    lv_image_t residulImage;
    residulImage.iNumAxes = inputImage->iNumAxes;
    residulImage.iPixelDepth = inputImage->iPixelDepth;
    residulImage.iFitsBitPix = inputImage->iFitsBitPix;
    residulImage.iTableDataType = inputImage->iTableDataType;
    residulImage.fPixelAngSize = inputImage->fPixelAngSize;
    residulImage.fDefaultVal = inputImage->fDefaultVal;
  */
    size_t imageSize = inputImage->pAxisSize[0]*inputImage->pAxisSize[1]; 
      for(int i=0; i<imageSize; ++i) {
          residulImage->pImage[i] = inputImage->pImage[i]-modelImage->pImage[i] ;
      }
       
    return residulImage;
    
  
}









/***************************
Function:	lv_read_img_file
Description:	Read a FITS file. Copied and modified from cookbook.c in cfitsio distribution
Arguments:	strFileName
Returns:	pointer to image stuct. NULL on failure.
****************************/
lv_image_t *lv_read_img_file(char	*strFileName, real_t fPixelRes) {

	fitsfile		*fptr=NULL;       /* pointer to the FITS file, defined in fitsio.h */
	lv_image_t		*pFitsImg = NULL;
	int				status=0,  nfound=0, anynull=0;
	lv_axissize_t	naxes[2], fpixel, nbuffer, npixels, ii, pixelType=0;
	char			strComment[1000];
	real_t			*pImg=NULL;

#define buffsize 1000
	real_t datamin, datamax, nullval, buffer[buffsize];

	TRACE_IN(lv_read_img_file);

	signal(SIGFPE, lv_fpesighandle);

	if (strFileName == NULL || strFileName[0] == '\0') {
		goto EXIT;
	}

    if ( fits_open_file(&fptr, strFileName, READONLY, &status) ) {
		sprintf(strMessage,"Failed to open file <%s>. exiting.",strFileName);
		TRACE(LOG_HIGH_PRI,strMessage);
		goto EXIT;
	}
	else {
		sprintf(strMessage,"Opened file <%s>",strFileName);
		TRACE(LOG_MED_PRI,strMessage);
	}

    /* read the NAXIS1 and NAXIS2 keyword to get image size */
    if ( fits_read_keys_lng(fptr, "NAXIS", 1,2, naxes, &nfound, &status) ) {
		TRACE(LOG_HIGH_PRI,"Failed to read key NAXIS. Exiting.");
		goto CLOSE_EXIT;
	}
	else {
		sprintf(strMessage,"Read axes: NAXIS1: %ld, NAXIS2: %ld.",(long) naxes[0],(long) naxes[1]);
		TRACE(LOG_MED_PRI,strMessage);
	}

	fits_read_key(fptr,TLONG, "BITPIX", &pixelType, strComment, &status);
	if (status != 0 ) {
		sprintf(strMessage,"Failed to read key BITPIX. status: %d. Exiting.",status);
		TRACE(LOG_HIGH_PRI,strMessage);
		goto CLOSE_EXIT;
	}
	else {
		sprintf(strMessage,"Read BITPIX. value: %ld",(long) pixelType);
		TRACE(LOG_MED_PRI,strMessage);
	}

    npixels  = naxes[0] * naxes[1];         /* number of pixels in the image */
    fpixel   = 1;
    nullval  = 0;                /* don't check for null values in the image */
    datamin  = FLT_MAX;
    datamax  = -FLT_MAX;

	/* create an image structure */
	pFitsImg = lv_create_image_struct(naxes[0], naxes[1], IMGTYPE, (fPixelRes == 0.0 ? g_PixelResn : fPixelRes));
	if (pFitsImg == NULL) {
		sprintf(strMessage,"lv_create_image_struct failed. exiting.");
		TRACE(LOG_HIGH_PRI,strMessage);
		goto CLOSE_EXIT;
	}
	pImg = pFitsImg->pImage;

	while (npixels > 0) {
		nbuffer = npixels;
		if (npixels > buffsize) {
			nbuffer = buffsize;
		}

		/* read as many pixels as will fit in buffer */
      /* Note that even though the FITS images contains unsigned integer */
      /* pixel values (or more accurately, signed integer pixels with    */
      /* a bias of 32768),  this routine is reading the values into a    */
      /* float array.   Cfitsio automatically performs the datatype      */
      /* conversion in cases like this.                                  */

/*
		if ( fits_read_img(fptr, TFLOAT, fpixel, nbuffer, &nullval, buffer, &anynull, &status) )
			printerror( status );
*/
		fits_read_img(fptr, pFitsImg->iTableDataType, fpixel, nbuffer, &nullval, buffer, &anynull, &status);
		if (status != 0) {
			sprintf(strMessage,"fits_read_img failed (status %d) after %ld pixels. exiting.", status, (long) fpixel-1);
			TRACE(LOG_HIGH_PRI,strMessage);
			lv_free_image_struct(pFitsImg);
			pFitsImg = NULL;
			goto CLOSE_EXIT;
		}
		memcpy(pImg, buffer, nbuffer * (size_t) pFitsImg->iPixelDepth);
		pImg += nbuffer;

		for (ii = 0; ii < nbuffer; ii++)  {
			if ( buffer[ii] < datamin )
				datamin = buffer[ii];

			if ( buffer[ii] > datamax )
				datamax = buffer[ii];
		}
		npixels -= nbuffer;    /* increment remaining number of pixels */
		fpixel  += nbuffer;    /* next pixel to be read in image */
	}

/*
    printf("\nMin and max image pixels =  %.0f, %.0f\n", datamin, datamax);
*/

CLOSE_EXIT:
    if ( fits_close_file(fptr, &status) ) {
		TRACE(LOG_HIGH_PRI,"Warning: fits_close_file failed. Continuing.");
	}

EXIT:
	signal(SIGFPE, SIG_DFL);
	TRACE_OUT;
    return pFitsImg;
}

static	void	lv_fpesighandle(int	sig) {
	fprintf(stderr,"Floating exception while loading FITS file. File contains illegal floating numbers.\n");
	exit(1);
}	


/***************************
Function:	lv_create_image_struct
Description:	create and initialise a 2D image structure
Arguments:
Returns:	pointer to the lv_image_t structure
****************************/
lv_image_t   *lv_create_image_struct(lv_axissize_t iXDimension, lv_axissize_t iYDimension, int iFitsPixelType,
            real_t fPixelAngSize) {
	lv_image_t	*pImg = NULL;
	
	pImg = (lv_image_t *) calloc(1,sizeof(lv_image_t));
	if (pImg == NULL) {
		return NULL;
	}
	pImg->iNumAxes = 2;
	pImg->iFitsBitPix = iFitsPixelType;
	pImg->fPixelAngSize = fPixelAngSize;
	pImg->iTableDataType = lv_decode_FITS_table_data_size(iFitsPixelType);
	pImg->iPixelDepth = lv_decode_FITS_pix_byte_size(iFitsPixelType);
	pImg->pAxisSize = (lv_axissize_t *) calloc(pImg->iNumAxes, sizeof(lv_axissize_t));
	if (pImg->pAxisSize == NULL) {
		free(pImg);
		return NULL;
	}
	pImg->pAxisSize[0] = iXDimension;
	pImg->pAxisSize[1] = iYDimension;
	pImg->pMask=NULL;
	pImg->bExternalMask = FALSE;
	pImg->ftInfo = NULL;
	pImg->pAxisOffset = calloc(pImg->iNumAxes,sizeof(float));
	if (pImg->pAxisOffset == NULL) {
		free(pImg->pAxisSize);
		free(pImg);
		return NULL;
	}
	pImg->pImage = (void *)calloc((iXDimension * iYDimension),pImg->iPixelDepth);
	if (pImg->pImage == NULL) {
		free(pImg->pAxisSize);
		free(pImg);
		return NULL;
	}

	return pImg;
}

/***************************
Function:	lv_duplicate_image
Description:	 Duplicate an image structure. Does not copy FT info.
Arguments:	pImg:	Image to be copied
Returns:	Pointer to new image struct, or NULL on failure.
****************************/
lv_image_t	*lv_duplicate_image(lv_image_t *pImg) {

	lv_image_t *pNewImg = NULL;
	size_t		i,iSize=0;

	if (pImg == NULL) {
		goto EXIT;
	}

	pNewImg = lv_create_image_struct(pImg->pAxisSize[0],pImg->pAxisSize[1],pImg->iFitsBitPix, pImg->fPixelAngSize);
	iSize = img_CalcImgSize(pImg);

	if (pNewImg == NULL) {
		goto EXIT;
	}

	lv_CopyRealArray(pImg->pImage,pNewImg->pImage,iSize);
	if (pImg->pMask != NULL) {
		pNewImg->pMask = (bool *) calloc(iSize,sizeof(bool));
		if (pNewImg->pMask == NULL) {
			lv_free_image_struct(pNewImg);
			pNewImg=NULL;
			goto EXIT;
		}
		for(i=0;i<iSize; i++) {
			pNewImg->pMask[i] = pImg->pMask[i];
		}
		pNewImg->bExternalMask = pImg->bExternalMask;
	}
	for (i=0; i< pImg->iNumAxes; i++) {
		pNewImg->pAxisOffset[i] = pImg->pAxisOffset[i];
	}
	pNewImg->fDefaultVal = pImg->fDefaultVal;

EXIT:
	return pNewImg;
}

/***************************
Function:
Description:
Arguments:
Returns:
****************************/
void	lv_free_image_struct(lv_image_t *pImagePtr) {

	if (pImagePtr != NULL) {
		free(pImagePtr->pAxisSize);
		free(pImagePtr->pAxisOffset);
		free(pImagePtr->pImage);
		if (pImagePtr->ftInfo != NULL) {
			free(pImagePtr->ftInfo);
		}
		if (pImagePtr->pMask != NULL) {
			free(pImagePtr->pMask);
		}
		free(pImagePtr);
	}
}


/*******************************
General Stuff
*******************************/
int		lv_decode_FITS_table_data_size(int iFitsPixelType) {
	int	iResult =0;

	switch (iFitsPixelType) {
		case BYTE_IMG :
			iResult = TBYTE;
			break;
		case SHORT_IMG:
			iResult = TSHORT;
			break;
		case LONG_IMG:
			iResult = TLONG;
			break;
		case FLOAT_IMG:
			iResult = TFLOAT;
			break;
		case DOUBLE_IMG:
			iResult = TDOUBLE;
			break;
	}
	return iResult;
}

/***************************
Function:	lv_ZeroRealArray
Description:	set all elements in a floating-point array to be zero.
Arguments:
Returns:
****************************/
void    lv_ZeroRealArray(register real_t *pArr, register lv_axissize_t iSize) {
	while(iSize-- > 0) {
		*(pArr++) = 0.0;
	}
}

/***************************
Function:	lv_CopyRealArray
Description:	copy iSize elements from pSrc to pDest
Arguments:
Returns:
****************************/
void	lv_CopyRealArray(register real_t *pSrc, register real_t *pDest, register lv_axissize_t iSize) {
	while(iSize-- > 0) {
		*(pDest++) = *(pSrc++);
	}
}

