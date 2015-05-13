/*******************
 * 	functions for calculating the critical lines
 * 	and caustics of lens models
 *	for the "lensview" software
 *	Copyright (C) Randall Wayth. 2006.
 *	$Id: lv_critline.c,v 1.4 2006/09/17 21:58:51 rwayth Exp $
*******************/
#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<string.h>
#include	"common.h"
#include	"log.h"
#include	"lv_common.h"
#include	"lv_geom.h"
#include	"lv_lens.h"
#include	"lv_critline.h"

/* private function prototypes */
static int	cl_CalcCritLinePixels(lv_image_t  *pImg, cl_pix_list *pList);
static int	cl_CritLinePixListAdd(cl_pix_list *pList,float x, float y);
static int	cl_CritLinePixListIsNeighbour(cl_pix_list *pList,float x, float y);
static int	cl_RefineCritLinePoints(lv_lensmodel_t *pLens, cl_pix_list *pList);

/***************************
Function:	cl_CalcCausticPoints
Description: Takes a list of points (actually pixels) which form the critical line and projects them
			in the source plane in the order then appear in the list. Forms a new set of pixels
			which are the points on the caustic.
			Assumes that the points in the crit line list form an orderly progression which
			traverses the critical line.
Arguments:	pLens: the current lens model.
			pCritlineList: list of pixels on the critical line. Only the x1,y1 values are used.
			pCausticList: output list of points on the caustic
Returns:	0: success. anything else: error
****************************/

int cl_CalcCausticPoints(lv_lensmodel_t *pLens, cl_pix_list *pCritlineList, cl_pix_list *pCausticList) {
	int iStatus=0;
	real_t	defly,deflx,mag=0;
	 cl_pix_list *pCurrent;

	TRACE_IN(cl_CalcCausticPoints);

	if (pCausticList==NULL || pCritlineList==NULL || pLens==NULL) {
		iStatus=1;
		LOG_ERR("input pointer NULL");
		goto EXIT;
	}

	pCurrent = NULL;

	while (pCritlineList != NULL) {
		/* if this isn't the (pre-allocated) first node in the list,
		 * then make a new node */
		if (pCurrent != NULL) {
			pCurrent->pNext = malloc(sizeof(cl_pix_list));
			if (pCurrent->pNext == NULL) {
				LOG_ERR("No malloc");
				iStatus = 1;
				goto EXIT;
			}
			pCurrent = pCurrent->pNext;
			pCurrent->pNext = NULL;
		}
		else {
			pCurrent = pCausticList;
		}

		/* must reset results */
		deflx = defly =0;

		iStatus = lm_CalcDeflection(pLens,pCritlineList->x1,pCritlineList->y1,&deflx,&defly,&mag);
		pCurrent->x1 = pCritlineList->x1 - deflx;
		pCurrent->y1 = pCritlineList->y1 - defly;
		sprintf(strMessage,"Crit line point (%g,%g) maps to src plane ( %g, %g )."
						,pCritlineList->x1,pCritlineList->y1,pCurrent->x1,pCurrent->y1);
		TRACE(LOG_MED_PRI,strMessage);
		pCritlineList = pCritlineList->pNext;
	}

EXIT:
	TRACE_OUT;
	return iStatus;
}

/***************************
Function:	cl_CalcCritLinePoints
Description:	Calculates an orderly succession of pixels on the critical line
				using a (pre-calculated) pixel inverse magnification map in the
				image plane.
Arguments:		pLens: current lens model. Must be the same as that used to generate
						the inv mag map.
				pImg: pixelised inverse magnification map for the image plane.
				pHead: placeholder for the head of the resulting list of critical line pixel values
Returns:	0: success. otherwise failure.
****************************/
int	cl_CalcCritLinePoints(lv_lensmodel_t *pLens, lv_image_t  *pImg, cl_pix_list *pHead) {
	int iStatus=0;
	int x,y,best_x=-1,best_y=-1,iHalfImgX,iHalfImgY;
	real_t	best_mag = 100;
	cl_pix_list *pTemp;

	TRACE_IN(cl_CalcCritLinePoints);

	/* find highest mag pixel. Start searching from there */
	for (x=0; x< pImg->pAxisSize[0]; x++) {
		for (y=0; y< pImg->pAxisSize[1]; y++) {
			if (pImg->pImage[x+y*pImg->pAxisSize[0]] < best_mag) {
				best_mag = pImg->pImage[x+y*pImg->pAxisSize[0]];
				best_x = x;
				best_y = y;
			}
		}
	}

	if (best_x == -1 || best_y == -1) {
		LOG_ERR("Failed to find any image plane pixels with inv mag < 100");
		iStatus=1;
		goto EXIT;
	}

	/* fill the placeholder head of the pixel list */
	pHead->pNext = NULL;
	pHead->x1 = best_x;
	pHead->x2 = best_x+1;
	pHead->y1 = best_y;
	pHead->y2 = best_y+1;

	/* create the actual list of pixels. This results in a rough critical line
	 * which is accurate to within a pixel width */
	iStatus = cl_CalcCritLinePixels(pImg,pHead);
	if (iStatus !=0) {
		LOG_ERR("cl_CalcCritLinePixels failed.");
		goto EXIT;
	}

	/* convert the pixel values to actual image plane x/y values
	 * in angular units */
	pTemp = pHead;
	iHalfImgX = pImg->pAxisSize[0]/2;
	iHalfImgY = pImg->pAxisSize[1]/2;
	while (pTemp != NULL) {
		pTemp->x1 = (pTemp->x1 - iHalfImgX)*pImg->fPixelAngSize;
		pTemp->y1 = (pTemp->y1 - iHalfImgY)*pImg->fPixelAngSize;
		pTemp->x2 = (pTemp->x2 - iHalfImgX)*pImg->fPixelAngSize;
		pTemp->y2 = (pTemp->y2 - iHalfImgY)*pImg->fPixelAngSize;
		pTemp = pTemp->pNext;
	}

	/* now refine the list by zooming in on each point more accurately */
	iStatus = cl_RefineCritLinePoints(pLens, pHead);

EXIT:
	TRACE_OUT;
	return iStatus;
}

/***************************
Function:	cl_RefineCritLinePoints
Description:	Takes a pre-calculated set of pixels which are the approximate
			critical line and refines the accuracy by subdividing the pixels
			until a fixel accuracy is reached.
Arguments:	pLens: the current lens model
			pList: List of pixels which lie on the critical line.
Returns:	0: success. otherwise failure.
****************************/
static int cl_RefineCritLinePoints(lv_lensmodel_t *pLens, cl_pix_list *pList) {
	int iStatus=0,i,j;
	bool bDone=FALSE;
	real_t imgx[3][3],imgy[3][3];
	real_t srcx[3][3],srcy[3][3];
	real_t x1,x2,y1,y2;
	real_t best_x1,best_x2,best_y1,best_y2,mag,best_mag,fArea1,fArea2;
	real_t deflx,defly,srcGridArea,srcTriAreaTot;
	lv_polygon SrcTri1,SrcTri2,SrcGrid;

	TRACE_IN(cl_RefineCritLinePoints);

	if (pLens== NULL || pList==NULL) {
		LOG_ERR("Null pointers");
		iStatus=1;
		goto EXIT;
	}

	/* initialise polygons used to calculate areas */
	SrcTri1.numVertices = 3;
	SrcTri2.numVertices = 3;
	SrcGrid.numVertices = 4;

	/* for each pixel in the list... */
	while (pList != NULL) {

		/* repeatedly subdivide pixel until we meet our uncertainty tolerance */
		bDone = FALSE;
		/* subdivide the pixel into 4 */
		x1 = pList->x1;
		x2 = pList->x2;
		y1 = pList->y1;
		y2 = pList->y2;
		while (!bDone) {

			/* fill temporary arrays with the x/y location of the 9 points in the
			 * subdivided pixel */
			for (i=0; i< 3; i++) {
				imgx[0][i]=x1;
				imgx[1][i]=(x1+x2)/2;
				imgx[2][i]=x2;
			}
			for (i=0; i<3; i++) {
				imgy[i][0]=y1;
				imgy[i][1]=(y1+y2)/2;
				imgy[i][2]=y2;
			}

			/* map each corner of each sub-divided pixel to the source plane */
			for(i=0; i< 3; i++) {
				for (j=0; j<3; j++) {
					deflx = defly =0;
					iStatus = lm_CalcDeflection(pLens,imgx[i][j],imgy[i][j],&deflx,&defly,&mag);
					srcx[i][j] = imgx[i][j] - deflx;
					srcy[i][j] = imgy[i][j] - defly;
					sprintf(strMessage,"Img (%g,%g) maps to src (%g,%g)",imgx[i][j],imgy[i][j],srcx[i][j],srcy[i][j]);
					TRACE(LOG_LOW_PRI,strMessage);
				}
			}

			/* form triangles in each subpixel and calculate the total area of the
			 * two projected triangles. The smallest area is the biggest mag
			 * so that pixel will contain the critical line */
			best_x1=0;
			best_y1=0;
			best_x2=0;
			best_y2=0;
			best_mag =100;
			srcTriAreaTot = 0;
			for(i=0; i< 2; i++) {
				for(j=0; j< 2 ; j++) {
					SrcTri1.vertex[1].x = srcx[i][j];
					SrcTri1.vertex[1].y = srcy[i][j];
					SrcTri1.vertex[2].x = srcx[i][j+1];
					SrcTri1.vertex[2].y = srcy[i][j+1];
					SrcTri1.vertex[3].x = srcx[i+1][j+1];
					SrcTri1.vertex[3].y = srcy[i+1][j+1];

					SrcTri2.vertex[1].x = srcx[i][j];
					SrcTri2.vertex[1].y = srcy[i][j];
					SrcTri2.vertex[2].x = srcx[i+1][j+1];
					SrcTri2.vertex[2].y = srcy[i+1][j+1];
					SrcTri2.vertex[3].x = srcx[i+1][j];
					SrcTri2.vertex[3].y = srcy[i+1][j];

					/* now calculate areas! */
					fArea1 = geomPolygonArea(&SrcTri1);
					fArea2 = geomPolygonArea(&SrcTri2);
					srcTriAreaTot += fArea1+fArea2;
					sprintf(strMessage,"Subregion (%d,%d) tri1, tri2 areas: %g,%g. Total: %g"
									,i,j,fArea1,fArea2,fArea1+fArea2);
					TRACE(LOG_LOW_PRI,strMessage);
					if (fArea1+fArea2 < best_mag) {
						best_x1 = imgx[i][j];
						best_x2 = imgx[i+1][j+1];
						best_y1 = imgy[i][j];
						best_y2 = imgy[j+1][j+1];
						best_mag = fArea1+fArea2;
					}
				}
			}

			/* sanity check: calculate the total area of the original pixel
			 * as projected into the source plane. If this pixel is actaully on
			 * the critical line,then it should be folded or twisted. Thus
			 * the area of this projected 4-sided polygon should not be the
			 * same as the sum of the areas from the projected triangles */
			SrcGrid.vertex[1].x = srcx[0][0];
			SrcGrid.vertex[1].y = srcy[0][0];
			SrcGrid.vertex[2].x = srcx[0][2];
			SrcGrid.vertex[2].y = srcy[0][2];
			SrcGrid.vertex[3].x = srcx[2][2];
			SrcGrid.vertex[3].y = srcy[2][2];
			SrcGrid.vertex[4].x = srcx[2][0];
			SrcGrid.vertex[4].y = srcy[2][0];

			srcGridArea = geomPolygonArea(&SrcGrid);

			/* check to see if this pixel is acutally folded */
			if (fabs((srcTriAreaTot-srcGridArea)/srcGridArea) < 1e-3) {
				sprintf(strMessage,"WARNING: Image grid does not appear to be folded. Grid area: %g. Total triangle area: %g",
								srcGridArea,srcTriAreaTot);
				LOG_ERR(strMessage);
			}

			/* have we got enough accuracy? */
			mag = sqrt(pow((best_x1-best_x2),2)+pow((best_y2-best_y1),2));
			if (mag < g_PixelResn/10) {
				bDone = TRUE;
			}
			sprintf(strMessage,"Updating (%g,%g),(%g,%g) to (%g,%g),(%g,%g)"
							,x1,y1,x2,y2,best_x1, best_y1,best_x2,best_y2);
			TRACE(LOG_LOW_PRI,strMessage);
			x1 = best_x1;
			x2 = best_x2;
			y1 = best_y1;
			y2 = best_y2;
		}

		/* put refined point back into pixel */
		sprintf(strMessage,"Refined crit line for pixel (%g,%g),(%g,%g) to (%g,%g),(%g,%g)"
						,pList->x1,pList->y1,pList->x2,pList->y2,x1,y1,x2,y2);
		TRACE(LOG_MED_PRI,strMessage);
		pList->x1 = x1;
		pList->x2 = x2;
		pList->y1 = y1;
		pList->y2 = y2;

		/* next pixel */
		pList = pList->pNext;
	}

EXIT:
	TRACE_OUT;
	return iStatus;
}

/***************************
Function:	cl_CalcCritLinePixels
Description:	Create an ordered list of pixels in the image plane
			which contain the critical line. Uses an already calculated
			pixel-by-pixel image of the inverse magnification.
			Also uses the first (and only) pixel in the list as the
			starting point.
Arguments:	pImg: image plane inverse magnification in each pixel
			pList: the output list of pixels
Returns:	0: success, otherwise failure
****************************/
static int	cl_CalcCritLinePixels(lv_image_t  *pImg, cl_pix_list *pList) {
	int iStatus=0,count=0,i,j,x,y,best_x,best_y,last_x=-1,last_y=-1;
	int	iMaxx,iMaxy;
	real_t	best_mag;
	bool	done=FALSE;

	TRACE_IN(cl_CalcCritLinePixels);

	iMaxx = pImg->pAxisSize[0] -1;
	iMaxy = pImg->pAxisSize[1] -1;

	/* init first pixel */
	x = pList->x1;
	y = pList->y1;

	sprintf(strMessage,"Beginning at pixel %d, %d",x,y);
	TRACE(LOG_MED_PRI,strMessage);

	while (!done) {

		sprintf(strMessage,"Searching for best mag around %d, %d. Last x:%d, y:%d",x,y,last_x,last_y);
		TRACE(LOG_LOW_PRI,strMessage);

		/* for each surrounding pixel, find the highest mag (lowest
		  pixel val) */
		best_mag=100;
		best_x = -1;
		best_y = -1;
		for (i=MAX(0,x-1); i<= MIN(iMaxx,x+1); i++) {
			for(j=MAX(0,y-1); j<= MIN(iMaxy,y+1); j++) {

				/* skip current pixel */
				if (i != x || j != y){
					/* if we are examining a pixel which is the same as
					 * the first pixel in the list, then we have come back
					 * to the start, so we are done */
					if(count > 10 && pList->x1==i && pList->y1==j) {
						done=TRUE;
					}
					else {
						int iIsNeighbour=0;

						/* check to see if this pixel is a neighbour of one already in the list.
						 * This stops us from going backwards. Must ignore the first and last
						 * items in the list since the last one is the previous pixel which must
						 * be a neighbour, and the first we check specially for stopping */
						iIsNeighbour = cl_CritLinePixListIsNeighbour(pList->pNext,i,j);
						sprintf(strMessage,"i: %d, j: %d, mag: %g, neighbour of exiting: %d"
										,i,j,pImg->pImage[i + j*pImg->pAxisSize[0]],iIsNeighbour);
						TRACE(LOG_LOW_PRI,strMessage);
						if (!(i == last_x && j == last_y) && pImg->pImage[i + j*pImg->pAxisSize[0]] < best_mag
										&& !iIsNeighbour) {
							best_x = i;
							best_y = j;
							best_mag = pImg->pImage[i + j*pImg->pAxisSize[0]];
							sprintf(strMessage,"best mag at %d, %d. mag: %g",i,j,best_mag);
							TRACE(LOG_LOW_PRI,strMessage);
						}
					}
				}
			}
		}
		if (last_x ==x && last_y ==y) {
			LOG_ERR("ERROR: got stuck following critical line. Stopping");
			done=TRUE;
		}
		if (!done) {
			sprintf(strMessage,"Adding pixel %d, %d to critline pixel list. Mag was %g. lastx: %d lasty: %d"
							,best_x,best_y,best_mag, last_x,last_y);
			TRACE(LOG_MED_PRI,strMessage);
			iStatus = cl_CritLinePixListAdd(pList,best_x,best_y);
			if (iStatus !=0 ) {
				LOG_ERR("CritLinePixListAdd failed");
				goto EXIT;
			}
		}
		last_x = x;
		last_y = y;
		x = best_x;
		y = best_y;
		count++;
	}
EXIT:
	TRACE_OUT;
	return iStatus;
}

/***************************
Function:	cl_CritLinePixListIsNeighbour
Description: check to see if a pixel position is within 1 pixel width
			of a pixel already in the list (including diagonal)
Arguments:	pList: current list of pixels
			x,y: position of pixel being examined
Returns:	0: not neighbour. 1: is neighbour
****************************/
static int	cl_CritLinePixListIsNeighbour(cl_pix_list *pList,float x, float y) {
	int count=0;
	
	if (pList == NULL) return 0;

	while(pList->pNext != NULL) {
		if (abs(pList->x1 - x) <=1 && abs(pList->y1 - y) <= 1) count++;
		pList = pList->pNext;
	}
	return count;
}

/***************************
Function:	cl_CritLinePixListAdd
Description:	Adds a new node to a list of pixels
Arguments:	pList: Head of the list
			x,y: pixel position
Returns:	0: success, otherwise failure
****************************/
static int	cl_CritLinePixListAdd(cl_pix_list *pList,float x, float y) {
	int iStatus=0;
	cl_pix_list *pTemp=NULL;

	TRACE_IN(cl_CritLinePixListAdd);

	while (pList->pNext != NULL) {
		pList = pList->pNext;
	}
	pTemp = malloc(sizeof(cl_pix_list));
	if (pTemp ==NULL) {
		LOG_ERR("No malloc");
		iStatus=1;
		goto EXIT;
	}
	pList->pNext = pTemp;
	pTemp->pNext = NULL;
	pTemp->x1 = x;
	pTemp->x2 = x+1;
	pTemp->y1 = y;
	pTemp->y2 = y+1;

EXIT:
	TRACE_OUT;
	return iStatus;
}


/***************************
Function:	cl_FreePixList
Description: Frees the malloced nodes of a pixel list
Arguments:	pList: start of the list to be freed
Returns:	0: success
****************************/
int	cl_FreePixList(cl_pix_list *pList) {
	int iStatus = 0;
	cl_pix_list *pTemp;

	while (pList != NULL) {
		pTemp = pList->pNext;
		free(pList);
		pList = pTemp;
	}

	return iStatus;
}

/***************************
Function:	cl_SavePixList
Description:	Write a pixel list to a text file as x, y values
			only the x1,y1 values of the pixels are written
Arguments:	pList: list of pixels
			strFilename: filename to be written to
Returns:	0: success, otherwise failure.
****************************/
int	cl_SavePixList(cl_pix_list *pList,char strFilename[]) {
	int iStatus = 0;
	FILE *fp;

	TRACE_IN(cl_SavePixList);
	if(pList == NULL || strFilename == NULL) {
		LOG_ERR("Null pointers");
		iStatus = 1;
		goto EXIT;
	}
	if ((fp=fopen(strFilename,"w"))==NULL) {
		sprintf(strMessage,"Cannot open output file <%s>",strFilename);
		LOG_ERR(strMessage);
		iStatus =1;
		goto EXIT;
	}

	while(pList != NULL) {
		fprintf(fp,"%g %g\n",pList->x1,pList->y1);
		pList = pList->pNext;
	}

	fclose(fp);
EXIT:
	TRACE_OUT;
	return iStatus;
}

