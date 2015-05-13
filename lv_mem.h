/****************
header files for the lv_mem.c module
part of the Lensview software package.
Copyright (C) 2006. Randall Wayth.
$Id: lv_mem.h,v 1.5 2008/10/27 01:08:04 rwayth Exp rwayth $
******************/

/*******************
Public funciton prototypes
*********************/
int mem_DoMEM(lv_lensmodel_t *pLens, lv_mapmatrix_t *pMap, lv_image_t *pSource, lv_image_t *pProjectedImage,
    lv_image_t *pMeasuredImage, int iMaxIterations, real_t *pFinalChiSqu, real_t *pFinalEntropy,
    lv_image_t **pBestSource, lv_image_t **pBestImage, lv_image_t *pPSF, lv_image_t *pNoise, real_t fSrcDefaultVal,
	real_t	*pFinalQ);


/*********************
 * global variables defined in this file
 * *******************/
extern	int g_bUseMultImgPixOnly;
