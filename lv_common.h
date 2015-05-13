#ifndef	LV_COMMON_H
#define	LV_COMMON_H

/*****************
lv_common.c/h: Funtions for commonly used/generic stuff in lensview
Copyright (C) 2004. Randall Wayth.
$Id: lv_common.h,v 1.20 2008/10/15 21:15:16 rwayth Exp rwayth $
******************/
#define	HST_WF_RESOLUTION		0.1	
#define	HST_PC_RESOLUTION		0.046
#define	HST_STIS_RESOLUTION		0.0507
#define	ARCSEC_IN_RADIAN		4.84813e-06
#define	MAX_FT_DIMENSION		4
#define	MAX_LENSCOMP_PARAMS		5

/*******************
Structure definitions
********************/

/* typedef placeholder for the axis size. FITS defines this as a long
 * which will cause problems if ints and longs aren't the same size */
typedef	long	lv_axissize_t;

/* type def for floating point numbers. If this is changed, then the
 * fits image type will need to be updated. Also, the typedef for points
 * in lv_geom.h will need to be updated */
 
/* this for double data */
/* 
*/
typedef	double	real_t;
#define	IMGTYPE	DOUBLE_IMG

/* this for float data
#define IMGTYPE	FLOAT_IMG
typedef float real_t;
 * */

/* floating point precision for map weights. This is where most of
 * the memory allocation in lensview goes, so floats will make
 * lensview use half as much memory as doubles. */
typedef float mapwt_t;

typedef struct _lv_image_t {
	short	iNumAxes;
	size_t	iPixelDepth;
	int		iFitsBitPix;
	int		iTableDataType;
	real_t	fPixelAngSize;
	real_t  fDefaultVal;        /* default source plane value */
	lv_axissize_t	*pAxisSize; /* X then Y ... */
	float	*pAxisOffset;		/* offset of the source plane from the optical axis */
	real_t	*pImage;	/* the actual data goes here */
	void	*ftInfo;
	bool	*pMask;		/* optional mask array. If mask value = TRUE, then use this pixel. If FALSE, ignore pixel*/
	bool	bExternalMask;
}	lv_image_t;

typedef	struct	_lv_lenscomponent_t {
	int		iType;
	int		iNumParams;
	real_t	fXoffset;
	real_t	fYoffset;
	real_t	fParameter[MAX_LENSCOMP_PARAMS];
	real_t	fParamfrom[MAX_LENSCOMP_PARAMS];
	real_t	fParamto[MAX_LENSCOMP_PARAMS];
	real_t	fParaminc[MAX_LENSCOMP_PARAMS];
	char	*strParamName[MAX_LENSCOMP_PARAMS];
	lv_image_t	*pDeflX;
	lv_image_t	*pDeflY;
}	lv_lenscomp;

typedef struct	_lv_lensmodel_t {
	int		iNumComponents;
	int		iNumParameters;
	lv_lenscomp	*pComponent;
}	lv_lensmodel_t;

typedef	struct	_lv_mapmatrix_t {
	lv_axissize_t	dimension[4];	/* Image x,y then Source x,y */
	mapwt_t			**array;
	lv_axissize_t	*min;
	lv_axissize_t	*max;
	real_t			*pSumImg;		/* the sum of the mapping matrix for a given image pixel */
	real_t			*pSumSrc;		/* the sum of the mapping matrix for a given source pixel weighted by
									   the sum of the image pixels mapping to it */
	bool			*pMultImgPix;	/* an array the size of the image. Each pixel represents whether the
									   image pixel is multiply imaged in the lens model or not */
}	lv_mapmatrix_t;

/*******************
Public Macro Definitions
*******************/
#define	lv_decode_FITS_pix_byte_size(a)	abs(a)/8
#define	NULL_PIXEL	0.0

/* NOTE: do not change types of global vars from float if they are affected by command line args */
extern	int		g_ShutdownFlag;
extern	float	g_FixedLambda;
extern	bool	*g_pSrcMask;
extern	float	g_DataMean;
extern	float	g_PixelResn;
extern	int     g_iDebugImgs;
extern	float	g_fTargetChiSqu;

/*******************
Public Function Prototypes
********************/

lv_mapmatrix_t	*lv_allocMappingMatrix(lv_axissize_t i,lv_axissize_t j,lv_axissize_t k,lv_axissize_t l);
int		lv_createMappingMatrix(lv_lensmodel_t *pLens, lv_mapmatrix_t *pMap, lv_image_t *pSourceImage, lv_image_t *pProjImage);
int		lv_freeMappingMatrix(lv_mapmatrix_t *pMap);
int		lv_zeroMapMatrix(lv_mapmatrix_t *pMap);

int		lv_projectRay(lv_lensmodel_t *pLens, real_t fImgX, real_t fImgY, real_t *pSrcX, real_t *pSrcY, real_t *pMag);
int		lv_projectImage(lv_lensmodel_t *lens, lv_image_t *source, lv_image_t *image);
int		lv_projectSourceThruMapMatrix(lv_image_t *pSource, lv_image_t *pImage, lv_mapmatrix_t *pMap);
int		lv_write_image_to_file(lv_image_t *pImagePtr, char  *strFileName, bool	bOverwriteFile);
lv_image_t	*lv_read_img_file(char   *strFileName, real_t fPixelRes);

lv_image_t	*lv_create_image_struct(lv_axissize_t iXDimension, lv_axissize_t iYDimension, int iFitsPixelType, real_t fPixelAngSize);
lv_image_t	*lv_duplicate_image(lv_image_t *pImage);
void    	lv_free_image_struct(lv_image_t *pImagePtr);

int		lv_decode_FITS_table_data_size(int iFitsPixelType);
void	lv_ZeroRealArray(register  real_t *pArr, lv_axissize_t iSize);
void	lv_CopyRealArray(real_t *pSrc, real_t *pDest, register lv_axissize_t iSize);

#endif	/* ifndef LV_COMMON_H */
