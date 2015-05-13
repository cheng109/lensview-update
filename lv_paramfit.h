/***************
lv_paramfit.h: header files for lv_paramfit.c
part of the Lensview software package.
Copyright (C) 2006. Randall Wayth.
$Id: lv_paramfit.h,v 1.2 2008/10/27 01:07:24 rwayth Exp rwayth $
*****************/
int lp_MinFinder(lv_lensmodel_t *pLens, lv_image_t *pSource, lv_image_t *pMappedImage, lv_image_t *pData,
                int iMaxIterations,real_t *fResChiSqu, real_t *fResEntropy, lv_image_t **pBestSrc, lv_image_t **pBestImg,
                lv_image_t *pPsfImg, lv_image_t *pNoiseImg, int bSearchCentre, float fSrcDefVal);

int lp_simplexWrap(lv_lensmodel_t *pLens, lv_mapmatrix_t *pMap, lv_image_t *pSource, lv_image_t *pProjectedImage,
            lv_image_t *pMeasuredImage, int iImgTraceFreq, real_t *pFinalChiSqu,
            real_t *pFinalEntropy, lv_image_t **pBestSource, lv_image_t **pBestImage, lv_image_t *pPSF,
            lv_image_t *pNoise, bool *pMask);
            
int lp_ParamSweep(lv_lensmodel_t *pLens, lv_image_t *pSourceOriginal, lv_image_t *pMappedImage, lv_image_t *pData,
                int iMaxIterations,real_t *pResChiSqu, real_t *pResEntropy, lv_image_t **pBestSrc,
                lv_image_t **pBestImg, lv_image_t *pPsfImg, lv_image_t *pNoiseImg, int bSearchCentre,
                float fSearchCentreRange, float fSearchCentreStep, float fSrcDefaultVal,
                int iMethod, int iDumpImgs);
