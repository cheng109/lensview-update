#ifndef LV_IMAGE_H
#define LV_IMAGE_H
/****************
header files for the lv_image.c module
part of the Lensview software package.
Copyright (C) 2006. Randall Wayth.
$Id: lv_image.h,v 1.17 2008/10/27 01:07:24 rwayth Exp rwayth $
******************/

/*******************
Public funciton prototypes
*********************/
int			img_AddArrays(real_t *pArr1, real_t *pArr2, real_t *pResult, size_t iSize);
lv_image_t	*img_convolveImages(lv_image_t *pImg, lv_image_t *pPSF, lv_image_t *pResult, bool *pMask);
int     img_ReverseTransposeImg(lv_image_t *pImg, lv_image_t *pResult);
real_t		img_CalcChiSquared(lv_image_t *pImg1, lv_image_t *pImg2, real_t fConstVariance, lv_image_t *pNoise, bool *pMask);
int			img_CalcImgEntropy(lv_image_t *pImg, real_t *pResult, real_t f_SkyVal);
int			img_DoConjGradient(lv_lensmodel_t *pLens, lv_mapmatrix_t *pMap, lv_image_t *pSource, lv_image_t *pProjectedImage,
				lv_image_t *pMeasuredImage, int iNumIterations, real_t *pFinalChiSqu,
				real_t *pFinalEntropy, lv_image_t **pBestSource, lv_image_t **pBestImage, lv_image_t *pPSF,
				lv_image_t *pNoise);
int		img_CalcEntropyDeriv(lv_image_t *pImg, real_t *pResult);
int		img_CalcChiSqDeriv(lv_image_t *pImg, lv_image_t *pdata, real_t *pResult, real_t ConstVariance, lv_image_t *pNoise, lv_image_t *pPSF);
int		img_CalcImgSize(lv_image_t *pImg);
int		img_CountNonZeroElements(register real_t *pArr, size_t iSize, real_t fThreshold);
int		img_CountActiveImgPixels(lv_image_t *pImg);
int		img_CountMultiplyImgdPixels(lv_mapmatrix_t *pMap);
real_t	img_MaxArrayVal(real_t *pArr, size_t iSize);
real_t	img_CalcArrayTotal(real_t *pArr, size_t iSize);
real_t	img_CalcVariance(real_t *pImg, lv_axissize_t iArrSize);
real_t	img_VecScalarProd(real_t *pArr1, real_t *pArr2, size_t iSize);
int		img_VecMultByScalar(real_t *pArr, real_t   fScalar, size_t iSize);
int		img_PadImage(lv_image_t *pImg, int *pDimension);
int     img_ConditionImage(lv_image_t *pImg, real_t fDefVal);
void	imgDumpMask(lv_image_t *pImg);
void    img_DumpImage(real_t *pImg, lv_axissize_t *pSizeChiArr, char *pName);
lv_image_t  *img_ConvolveImgsWithFFT(lv_image_t *pImg, lv_image_t *pPSF, lv_image_t *pResultImage);
lv_image_t	*img_CalcSrcPlaneMag(lv_mapmatrix_t *pMap, real_t fSrcPixAngSize, real_t fImgPixAngSize);
lv_image_t  *img_CalcImgPlaneInvMag(lv_mapmatrix_t *pMap, real_t fSrcPixAngSize, real_t fImgPixAngSize);
int     img_ReverseProject(lv_mapmatrix_t *pMap, lv_image_t *pChiDeriv, real_t *pSourceDelta);
int		img_FindMultiplyImgdPixels(lv_mapmatrix_t *pMap);

/*******************
Public global variables
*********************/
extern	real_t  g_imgImgVariance;
extern	bool	g_imgCalcVarianceFlag;

#endif	/* LV_IMAGE_H */
