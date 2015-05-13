/******************************
lv_image.c:	image/array manipulation and lensing-specific optimisation
	programs for the "lensview" software
Copyright (C) 2006. Randall Wayth
$Id: lv_image.c,v 1.37 2008/10/31 20:43:10 rwayth Exp rwayth $
*******************************/
#include	<unistd.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>
#include	<time.h>
#include	<limits.h>
#include	<float.h>
#include	"log.h"
#include	"common.h"
#include	"lv_common.h"
#include	"lv_image.h"
#include	"fitsio.h"
/* single-precision FFT
#include	"srfftw.h"
*/
/* double-precision FFT
#include	"drfftw_threads.h"
*/
#include	"rfftw.h"

#define		SMALL_VAL	1E-7
#define		MAX_CHI_TOLERANCE	10
#define		MIN_CHI_IMPROV		0.01
#define		EPS	3.0e-8
#define		ZEPS	1e-20
#define		GOLD_RATIO		1.618034
#define		GOLD_RATIO_INV	0.3819660
#define		MAG_LIMIT		100.0
#define		SHFT(a,b,c,d)	(a)=(b);(b)=(c);(c)=(d);
#define		SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define		MIN_TOL	1e-4
#define		ITMAX	100
#define		UPDATE_LAMBDA	-2
#define		LAMDA_START	5.0
#define		LAMBDA_TOL	0.1

/*************************
Public Global Variables
*************************/
real_t	g_imgImgVariance = 1.0;
bool	g_imgCalcVarianceFlag = TRUE;

/*************************
Private Global Variables
*************************/
static	real_t	g_Lambda = LAMDA_START;
static	int		g_iSizeImg=0;
static	int		g_iSizeSrc=0;
static	real_t	*g_pCommonImg=NULL, *g_pCommonSrc=NULL, *g_pXiImg=NULL, *g_pXiSrc=NULL;
static	lv_image_t	*g_pData=NULL, *g_pVariance=NULL, *g_pTempImg=NULL, *g_pTempSrc=NULL;
static	lv_image_t	*g_pPSF=NULL, *g_pRevPSF=NULL;
static	lv_mapmatrix_t	*g_pMapMatrix=NULL;

/*************************
Private Functions
*************************/
static	int		img_Calc1dJ(real_t x, real_t *pResultChiSqu, real_t *pResultEnt);
static	int		img_CalcJ(lv_image_t *pImg, lv_image_t *pSrc, real_t *pResChiSqu, real_t *pResEntropy, lv_image_t *pNoise, lv_image_t *pPSF);
static	int		img_CalcGradJ(lv_image_t *pImg, lv_image_t *pSrc, real_t *pGradSrc);
static	int		img_Calc1dGradJ(real_t x, real_t *pResult);
static	int		img_Bracket1dMin(real_t *a, real_t *b, real_t *c, real_t *fa, real_t *fb, real_t *fc);
static	int		img_DoBrentDeriv(real_t brack_a, real_t brack_b, real_t brack_c, real_t fTolerance, real_t *pXmin, real_t *pJmin);
static	int		img_DoBrent(real_t brack_a, real_t brack_b, real_t brack_c, real_t fTolerance, real_t *pXmin, real_t *pChiSqu, real_t *pEntropy);
static	int		img_LineMin(lv_image_t *pImg, lv_image_t *pSrc, real_t *pdelSrc, real_t *pResChiSq, real_t *pResEnt);
int img_Minimise(lv_lensmodel_t *pLens, lv_mapmatrix_t *pMap, lv_image_t *pSource, lv_image_t *pImg,
    lv_image_t *pData, int iMaxIterations, real_t fTolerance, real_t *pFinalChiSq, real_t *pFinalEnt,
    lv_image_t *pPSF, lv_image_t *pNoise, int *iIter);

/***************************
chi-sq minimisation

Many routines here based on optimisation routines from Numerical Recipes in C 2nd Edition chapter 10.
They have been customised for grav lensing and made error-friendly.

****************************/

/***************************
Function:	img_CalcJ
Description: Calculate J = Entropy - lambda*chisquare for image and source
			does not assume that the image has been constructed. Constructs the
			image by projecting and convolving with PSF and puts result in pImg
Arguments:	pImg, pSrc. Image and source. pResult resulting value
Returns:	0: success. Non-zero: failure
****************************/
static	int	img_CalcJ(lv_image_t *pImg, lv_image_t *pSrc, real_t *pResultChiSqu, real_t *pResEntropy, lv_image_t *pNoise, lv_image_t *pPSF){
	int	iRes=0;

	TRACE_IN(img_CalcJ);

	/* condition the source so there are no negative values */
	if (g_Lambda !=0) {
		img_ConditionImage(pSrc,SMALL_VAL);
	}

	/* project the new source through the map matrix */
	lv_projectSourceThruMapMatrix(pSrc, pImg, g_pMapMatrix);

	/* convolve with the PSF */
	if (pPSF != NULL) {
		sprintf(strMessage,"Convolving mapped image with PSF.");
		TRACE(LOG_LOW_PRI,strMessage);

		/* NOTE: the rotated psf is used here so that the convolution
		 * simulated the true effects of a psf. See note in function
		 * definition */

		g_pTempImg = img_ConvolveImgsWithFFT(pImg,pPSF,g_pTempImg);
		lv_CopyRealArray(g_pTempImg->pImage,pImg->pImage,g_iSizeImg);
	}

	/* calc chi-squared. 1st arg must be data. 2nd must be model */
	*pResultChiSqu = img_CalcChiSquared(g_pData,pImg, (pNoise==NULL ? g_imgImgVariance : 0), pNoise, pImg->pMask);

	/* calc Entropy */
	/* for the special case when lambda=0, we allow negative source pixels so can't calc entropy */
	if (g_Lambda !=0) {
		img_CalcImgEntropy(pSrc,pResEntropy, g_DataMean);
	}
	else {
		*pResEntropy=0;
	}

	sprintf(strMessage,"chi-sq is: %g, entropy is: %g",*pResultChiSqu,*pResEntropy);
	TRACE(LOG_LOW_PRI,strMessage);


	TRACE_OUT;
	return iRes;
}

/***************************
Function:	img_Calc1dJ
Description:	Calculates J using funciton above. The source is constucted using
				a scalar multiplication factor for a "delta" to the source. The
				new source is construcuted as: s_new = s + x*s_delta
				uses global variables g_pCommonSrc, g_pXiSrc which are the original
				source and source delta respectively.
Arguments:	x: scale factor, pResChiSqu: resulting value of chi-squ, pResEntropy: resulting source entropy
Returns:	0: success. Non-zero: failure
****************************/
int	img_Calc1dJ(real_t x, real_t *pResChiSqu, real_t *pResEntropy){
	static	lv_image_t	*pTempImg = NULL, *pTempSrc = NULL;

	int iRes=0,i;

	TRACE_IN(img_Calc1dJ);

	if (pTempImg == NULL) {
		pTempImg = lv_duplicate_image(g_pTempImg);	/* note: must use g_pTempImg here because it has a mask array */
		pTempSrc = lv_duplicate_image(g_pTempSrc);
		if (pTempImg == NULL || pTempSrc == NULL) {
			LOG_ERR("lv_duplicate_image failed. No malloc.");
			iRes = 1;
			goto EXIT;
		}
	}

	/* construct a temp img/src given the current scale factor on the gradient x */
	for(i=0; i<g_iSizeSrc; i++) {
		pTempSrc->pImage[i] = g_pCommonSrc[i] + x*g_pXiSrc[i];
	}

	iRes = img_CalcJ(pTempImg,pTempSrc,pResChiSqu,pResEntropy,g_pVariance,g_pPSF);

	sprintf(strMessage,"chi-squ is %g, Ent is %g for gradient scaling of %g",*pResChiSqu,*pResEntropy,x);
	TRACE(LOG_LOW_PRI,strMessage);

EXIT:
	TRACE_OUT;
	return iRes;
}


/***************************
Function:	img_CalcGradJ
Description:	Calculates a 2-part "gradient" vector for the source/image combination
			Reverse-projects the lambda-scaled gradient from the image plane and adds
			to the gradient from the source to create a single gradient vector which is
			a weighted combination of the gradients from source and image.
			Uses global variables for data, psf, image sizes etc.
Arguments:	pImg, pSrc: Image and source. pGradSrc: vector for resulting combined gradient
Returns:	0: success. Non-zero: failure.
****************************/
static	int	img_CalcGradJ(lv_image_t *pImg, lv_image_t *pSrc, real_t *pGradSrc) {
	int	iRes =0,i=0;

	TRACE_IN(img_CalcGradJ);

	/* calc the gradient from the image. Put the result in g_pTempImg. Result is already PSF convolved */
	iRes = img_CalcChiSqDeriv(g_pData,pImg,g_pTempImg->pImage,(g_imgCalcVarianceFlag == TRUE ? 0.0 : g_imgImgVariance), g_pVariance, g_pRevPSF);
	if (iRes != 0) {
		goto EXIT;
	}
	if (g_iDebugImgs)  {
		img_DumpImage(g_pTempImg->pImage,g_pTempImg->pAxisSize,"gradC");
	}


	/* calc the gradient from the source. Put the result in pGradSrc */
	/* in the special case of lambda=0, we allow negative source pixels so can't calc
		gradient (not that it matters because we mult by 0 later) */
	if (g_Lambda !=0) {
	    pSrc->fDefaultVal=g_DataMean;
		iRes = img_CalcEntropyDeriv(pSrc,pGradSrc);
		if (iRes != 0) {
			goto EXIT;
		}
	}
	else {
		lv_ZeroRealArray(pGradSrc,g_iSizeSrc);
	}

	/* reverse project the PSF-convolved gradient vector back to src */
	/* put the result in g_pTempSrc */
	lv_ZeroRealArray(g_pTempSrc->pImage, g_iSizeSrc);

	iRes = img_ReverseProject(g_pMapMatrix, g_pTempImg, g_pTempSrc->pImage);
	if (iRes != 0) {
		goto EXIT;
	}

	if (g_iDebugImgs)  {
		img_DumpImage(g_pTempSrc->pImage,g_pTempSrc->pAxisSize,"revprojGradC");
	}

	/* scale the entropy gradient by Lambda */
	img_VecMultByScalar(pGradSrc,g_Lambda,g_iSizeSrc);

	/* mult the chi gradient by -1 because that is what is needed later */
	for (i=0; i< g_iSizeSrc; i++) g_pTempSrc->pImage[i] = -(g_pTempSrc->pImage[i]);

	/* add this reverse projected image to the source which already contains
	 * the gradient from the entropy */
	iRes = img_AddArrays(g_pTempSrc->pImage, pGradSrc, pGradSrc, g_iSizeSrc);

EXIT:
	TRACE_OUT;
	return iRes;
}


/***************************
Function:
Description:
Arguments:
Returns:
****************************/
/************
	NOT USED. PROBABLY DOES NOT WORK
************/
static	int	img_Calc1dGradJ(real_t x, real_t *pResult) {

	static	lv_image_t	*xTempSrc=NULL;
	static	real_t	*dfSrc=NULL;

	register	int	i;
	register	real_t	df1 =0.0;
	int			iRes=0;

	TRACE_IN(img_Calc1dGradJ);

	if (xTempSrc==NULL) {
		xTempSrc = lv_duplicate_image(g_pTempSrc);
		dfSrc = malloc(sizeof(real_t)*g_iSizeSrc);
		if (dfSrc == NULL || xTempSrc == NULL) {
			LOG_ERR("image duplication/alloc failed. No malloc.");
			iRes=1;
			goto EXIT;
		}
	}
	if (iRes != 0) {
		goto EXIT;
	}

	for(i=0; i< g_iSizeSrc; i++) {
		df1 += dfSrc[i]*xTempSrc->pImage[i];
	}
	
	*pResult = df1;
EXIT:
	TRACE_OUT;
	return iRes;
}


/***************************
Function:	img_Bracket1dMin	 (based on NumRec function "mnbrak")
Description:	Find a region which bounds a minimum for the funciton J.
Arguments:	a,b,c: initial guess/result for bounding regions. fa,fb,fc: values of J
			at the corresponding points for a,b,c.
Returns:	0: success. non-zero: failure.
****************************/
static	int		img_Bracket1dMin(real_t *a, real_t *b, real_t *c, real_t *fa, real_t *fb, real_t *fc) {
	real_t	ulim,u,r,q,fu,dummy,qminr,tempchi,tempent;
	int	iRes=0,count=0;

	TRACE_IN(img_Bracket1dMin);

	iRes = img_Calc1dJ(*a,&tempchi,&tempent);
	if (iRes != 0) {
		goto EXIT;
	}
	*fa = tempchi - g_Lambda*tempent;

	/* if non-zero delta-source is larger than the starting
	 * point, then we probably have a too big guess for 
	 * the initial b, so try to reduce it */
	do {
		*b /= 10;
		iRes = img_Calc1dJ(*b,&tempchi,&tempent);
		if (iRes != 0) {
			goto EXIT;
		}
		*fb = tempchi - g_Lambda*tempent;
		sprintf(strMessage,"Looking for J < %g. fb is %g for b %g",*fa,*fb,*b);
		TRACE(LOG_MED_PRI,strMessage);
	} while ((*fb > *fa) && ++count < 10);
	if (count >= 10) {
		LOG_ERR("Can't find small enough initial b for bracket.");
		iRes=1;
		goto EXIT;
	}

	sprintf(strMessage,"Beginning with bounds (%g,%g) with J values (%g,%g)",*a, *b, *fa,*fb);
	TRACE(LOG_MED_PRI,strMessage);

	/* make sure we are going downhill from a to b */
	if (*fb > *fa) {
		SHFT(dummy,*a,*b,dummy)
		SHFT(dummy,*fa,*fb,dummy)
	}

	/* first guess for c */
	*c = *b + GOLD_RATIO*(*b - *a);

	iRes = img_Calc1dJ(*c,&tempchi,&tempent);
	if (iRes != 0) {
		goto EXIT;
	}
	*fc = tempchi - g_Lambda*tempent;

	while (*fb > *fc) {

		sprintf(strMessage,"a: %g, b: %g, c: %g, fa: %g, fb: %g, fc: %g",*a,*b,*c,*fa,*fb,*fc);
		TRACE(LOG_MED_PRI,strMessage);

		r = (*b - *a)*(*fb - *fc);
		q = (*b - *c)*(*fb - *fa);
		qminr= q-r;

		/* obvious, isn't it? */
		u = (*b)-((*b -*c)*q-(*b-*a)*r)/(2.0*SIGN(MAX(fabs(qminr),1e-20),qminr));
		ulim = (*b) + MAG_LIMIT*(*c-*b);

		/* try parabolic fit between b and c */
		if((*b-u)*(u-*c) > 0.0) {
			iRes = img_Calc1dJ(u,&tempchi,&tempent);
			if (iRes != 0) {
				goto EXIT;
			}
			fu = tempchi - g_Lambda*tempent;
			if(fu < *fc) {
				/* min between b and c */
				*a = *b;
				*b = u;
				*fa = *fb;
				*fb = fu;
				goto EXIT;
			}
			else if (fu > *fb) {
				/* min between a and u */
				*c = u;
				*fc = fu;
				goto EXIT;
			}

			/* parabolic fit didn't work. Use default golden ratio */
			u = (*c)+GOLD_RATIO*(*c - *b);
			iRes = img_Calc1dJ(u,&tempchi,&tempent);
			if (iRes != 0) {
				goto EXIT;
			}
			fu = tempchi - g_Lambda*tempent;
		}
		/* try between c and limit */
		else if ((*c -u)*(u - ulim) > 0.0) {
			iRes = img_Calc1dJ(u,&tempchi,&tempent);
			if (iRes != 0) {
				goto EXIT;
			}
			fu = tempchi - g_Lambda*tempent;
			if (fu < *c) {
				real_t	fTemp;
				SHFT(*b,*c,u,*c+GOLD_RATIO*(*c - *b))
				iRes = img_Calc1dJ(u,&tempchi,&tempent);
				if (iRes != 0) {
					goto EXIT;
				}
				fTemp = tempchi - g_Lambda*tempent;
				SHFT(*fb,*fc,fu,fTemp)
			}
		}
		/* limit u to max allowed value */
		else if((u-ulim)*(ulim-*c) >= 0.0) {
			u = ulim;
			iRes = img_Calc1dJ(u,&tempchi,&tempent);
			if (iRes != 0) {
				goto EXIT;
			}
			fu = tempchi - g_Lambda*tempent;
		}
		/* use golden method */
		else {
			u = (*c)+GOLD_RATIO*(*c - *b);
			iRes = img_Calc1dJ(u,&tempchi,&tempent);
			if (iRes != 0) {
				goto EXIT;
			}
			fu = tempchi - g_Lambda*tempent;
		}
		/* throw out oldest point and continue */
		SHFT(*a,*b,*c,u)
		SHFT(*fa,*fb,*fc,fu)
	}

EXIT:
	sprintf(strMessage,"status %d. a: %g, b: %g, c: %g, fa: %g, fb: %g, fc: %g", iRes,*a,*b,*c,*fa,*fb,*fc);
	TRACE(LOG_LOW_PRI,strMessage);

	TRACE_OUT;
	return iRes;
}


/***************************
Function:	img_DoBrent	(based on NumRec function "brent")
Description:	Use Brent's method to find a minimum in the 1-dimentional function J which is bounded
				by a and c up to a given tolerance. a and c are really the scale factors for the
				gradient vector which is used in img_Calc1dJ
Arguments:	a,c: bounds. b: initial guess. fTolerance: min relative improvement in J before done.
			pXmin: pXmin scale factor value which produces minimum. pJmin: corresponding J value at min.
Returns:	0: success. non-zero: failure
****************************/
static	int	img_DoBrent(real_t brack_a, real_t brack_b, real_t brack_c, real_t fTolerance, real_t *pXmin, real_t *pChiSqu, real_t *pEnt) {
	real_t	a,b,d=0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm,e=0.0,fTempChi=0,fTempEnt=0;
	int	iRes=0,i;

	TRACE_IN(img_DoBrent);

	sprintf(strMessage,"Beginning with ax,bx,cx: %g,%g,%g. Tol: %g",brack_a,brack_b,brack_c,fTolerance);
	TRACE(LOG_MED_PRI,strMessage);

	a = MIN(brack_a,brack_c);
	b = MAX(brack_a,brack_c);

	/* init and calculate inital value of J */
	x = w = v = brack_b;
	iRes = img_Calc1dJ(x,&fTempChi,&fTempEnt);
	if (iRes != 0) {
		goto EXIT;
	}
	fw = fv = fx = (fTempChi - g_Lambda*fTempEnt);

	for (i=0; i< ITMAX; i++) {

		/* test for doneness */
		xm = 0.5*(a+b);
		tol1 = fTolerance * fabs(x)+ZEPS;
		tol2 = 2.0*tol1;
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*pXmin = x;
			*pChiSqu = fTempChi;
			*pEnt = fTempEnt;
			sprintf(strMessage,"Done. (iter=%d). x: %g, xm: %g, a: %g, b: %g, chi: %g, ent: %g, J: %g",i,x,xm,a,b,fTempChi,fTempEnt,fTempChi - g_Lambda*fTempEnt);
			TRACE(LOG_MED_PRI,strMessage);
			goto EXIT;
		}

		if (fabs(e) > tol1) {
			r = (x-w)*(fx-fv);
			q = (x-v)*(fx-fw);
			p = (x-v)*q - (x-w)*r;
			q = 2.0 * (q-r);
			if (q > 0.0) {
				p = -p;
			}
			q = fabs(q);
			etemp = e;
			e=d;
			if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x)) {
				e = x >= xm ? a-x : b-x;
				d = GOLD_RATIO_INV*e;
			}
			else {
				d = p/q;
				u = x+d;
				if (u-a < tol2 || b-u < tol2) {
					d = SIGN(tol1,xm-x);
				}
			}
		}
		else {
			e = x >= xm ? a-x : b-x;
			d = GOLD_RATIO_INV*e;
		}
		u = fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d);
		iRes = img_Calc1dJ(u,&fTempChi,&fTempEnt);
		if (iRes != 0) {
			goto EXIT;
		}
		fu = fTempChi- g_Lambda*fTempEnt;
		if (fu <= fx) {
			if (u >= x) {
				a=x;
			}
			else {
				b=x;
			}
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		}
		else {
			if (u < x) {
				a=u;
			}
			else {
				b=u;
			}
			if (fu <= fw || w == x) {
				v=w; w=u; fv=fw; fw=fu;
			}
			else if(fu <= fv || v==x || v==w) {
				v=u; fv=fu;
			}
		}
	}
	LOG_ERR("Too many iterations");
	iRes = 1;

EXIT:
	TRACE_OUT;
	return iRes;
}


/***************************
Function:
Description:	**** NOT USED PROBABLY DOESN'T WORK.*****
Arguments:
Returns:
****************************/
static	int	img_DoBrentDeriv(real_t brack_a, real_t brack_b, real_t brack_c, real_t fTolerance, real_t *pXmin, real_t *pJmin) {
	int	i,ok1,ok2, iRes=0;
	real_t	a,b,d,d1,d2,du,dv,dw,dx,e=0.0,fTempChi,fTempEnt;
	real_t	fx,fu,fv,fw,old_e=0,tol1,tol2,u,u1,u2,v,w,x,xm;

	TRACE_IN(img_DoBrentDeriv);

	a = MIN(brack_a,brack_c);
	b = MAX(brack_a,brack_c);
	x = w = v = brack_b;

	iRes = img_Calc1dJ(x,&fTempChi,&fTempEnt);
	if (iRes != 0) {
		goto EXIT;
	}
	fw = fv = fx = fTempChi;
	iRes = img_Calc1dGradJ(x,&fTempChi);
	if (iRes != 0) {
		goto EXIT;
	}
	dw = dv = dx = fTempChi;

	for(i=0; i<ITMAX; i++) {
		xm = 0.5*(a+b);
		tol1 = fTolerance*fabs(x) + EPS;
		tol2 = 2.0 * tol1;

		sprintf(strMessage,"Beginning iteration %d. a: %g, b: %g, x: %g, w: %g, v: %g, fx: %g, fw: %g, fv: %g",i,a,b,x,w,v,fx,fw,fv);
		TRACE(LOG_MED_PRI,strMessage);

		if(fabs(x-xm) <= (tol2 - 0.5*(b-a))) {
			*pXmin = x;
			*pJmin = fx;
			goto EXIT;
		}

		/* make parabolic fits */
		if(fabs(e) > tol1) {
			d2 = d1 = 2.0*(b-a);
			if (dw != dx) {
				d1 = (w-x)*dx/(dx-dw);
			}
			if (dv != dx) {
				d2 = (v-x)*dx/(dx-dv);
			}
			u1 = x+d1;
			u2 = x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			old_e = e;
			e = d;

			/* check to see if either trial is OK. Take the smallest if both */
			if (ok1 || ok2) {
				if (ok1 && ok2) {
					d = (fabs(d1) < fabs(d2) ? d1 : d2);
				}
				else if (ok1){
					d = d1;
				}
				else {
					d = d2;
				}

				/* check that the step size is small enough */
				if (fabs(d) <= fabs(0.5*old_e)) {
					u = x+d;
					if (u-a < tol2 || b-u < tol2) {
						d = SIGN(tol1,xm-x);
					}
				}
				/* otherwise just bisect the points */
				else {
					d = 0.5 * (e=(dx >= 0.0 ? a-x : b-x));
				}
			}
			else {
				d = 0.5 * (e=(dx >= 0.0 ? a-x : b-x));
			}
		}
		else {
			d = 0.5 * (e=(dx >= 0.0 ? a-x : b-x));
		}

		/* d contains the best offset for x? */
		if (fabs(d) >= tol1) {
			u=x+d;
			iRes = img_Calc1dJ(u,&fTempChi,&fTempEnt);
			fu = fTempChi- g_Lambda*fTempEnt;
		}
		else {
			u = x+SIGN(tol1,d);
			iRes = img_Calc1dJ(u,&fTempChi,&fTempEnt);
			fu = fTempChi- g_Lambda*fTempEnt;
			if (fu > fx) {
				*pXmin =x;
				*pJmin = fx;
				goto EXIT;
			}
		}

		/* calculate new gradient and update housekeeping variables */
		iRes = img_Calc1dGradJ(u,&fTempChi);
		du = fTempChi- g_Lambda*fTempEnt;
		if(fu <= fx) {
			if (u >= x) {
				a=x;
			}
			else {
				b=x;
			}
			v=w; fv = fw; dv = dw;
			w=x; fw = fx; dw = dx;
			x=u; fx = fu; dx = du;
		}
		else {
			if (u < x) {
				a=u;
			}
			else {
				b=u;
			}
			if (fu <= fw || w==x) {
				v=w; fv = fw; dv = dw;
				w=u; fw = fu; dw = du;
			}
			else if (fu < fv || v == x || v == w) {
				v=u; fv = fu; dv = du;
			}
		}
	}
	LOG_ERR("Too many iterations")
	iRes =1;
EXIT:
	TRACE_OUT;
	return iRes;
}


/***************************
Function:	img_LineMin	(based on NumRec "dlinmin")
Description:	Do 1-d line minimisation on the function J in the search direction pdelSrc.
Arguments:	pImg,pSrc: Initial image source and destination for results.
			pdelSrc: gradient vector for source plane which will be searched along.
			pResMinJ: resulting minimum value of J in direction of gradient
Returns:	0: success. non-zero: failure
****************************/
static	int	img_LineMin(lv_image_t *pImg, lv_image_t *pSrc, real_t *pdelSrc, real_t *pResChiSqu, real_t *pResEntropy) {
	int i,iRes=0;
	lv_image_t	*pImgDelta = NULL;
	real_t	x, xMin =1.0,fx=0,fa=0,fb=0,b=0,a=0, fTemp;

	TRACE_IN(img_LineMin);

	/* set up the common data shared by this and img_Calc1dJ */
	/* these should *NOT* be changed. Updates are made in Minimising function */
	/* changes to these vals for best fit attempts should be done in new
	 * working space */
	g_pCommonImg = pImg->pImage;
	g_pCommonSrc = pSrc->pImage;
	g_pXiSrc = pdelSrc;
	pImgDelta = lv_duplicate_image(pImg);
	if(pImgDelta == NULL) {
		LOG_ERR("image duplicate failed. no malloc.");
		iRes = -1;
		goto EXIT;
	}
	lv_ZeroRealArray(pImgDelta->pImage,g_iSizeImg);

	/* project the source gradient through the map matrix */
	/* this will give us a "static" image which can be added to the
	 * current best img for calculating J with scaled src/projected src */
	lv_CopyRealArray(pdelSrc,g_pTempSrc->pImage,g_iSizeSrc);
	lv_projectSourceThruMapMatrix(g_pTempSrc, pImgDelta, g_pMapMatrix);

	/* convolve with the PSF */
	if (g_pPSF != NULL) {
		sprintf(strMessage,"Convolving mapped image with PSF.");
		TRACE(LOG_MED_PRI,strMessage);
		g_pTempImg = img_ConvolveImgsWithFFT(pImgDelta,g_pPSF,g_pTempImg);
		lv_CopyRealArray(g_pTempImg->pImage,pImgDelta->pImage,g_iSizeImg);
	}
	g_pXiImg = pImgDelta->pImage;

	/* initial guesses for bounds */
	a = 0.0;
	x = 1.0;
	fTemp = img_MaxArrayVal(pdelSrc,g_iSizeSrc);
	if (fabs(fTemp) > 1.0) {
		x = 0.1/fTemp;
	}
	else {
		x = 0.1*fTemp;
	}

	/* find a bounding region which must have a minimum */
	/* b is set inside function */
	iRes = img_Bracket1dMin(&a,&x,&b,&fa,&fx,&fb);
	if (iRes != 0) {
		LOG_ERR("img_Bracket1dMin failed. exiting.");
		goto EXIT;
	}

	sprintf(strMessage,"Found min bounds. a: %g, x: %g, c: %g, fa: %g, fx: %g, fc: %g",a,x,b,fa,fx,fb);
	TRACE(LOG_MED_PRI,strMessage);

	if (x < 0) {
		sprintf(strMessage,"img_Bracket1dMin returned a negative scale factor %g.",x);
		TRACE(LOG_MED_PRI,strMessage);
		iRes =UPDATE_LAMBDA;
		goto EXIT;
	}

	/* now find the actual minimum */
	iRes = img_DoBrent(a,x,b,1e-2,&xMin,pResChiSqu,pResEntropy);
	if (iRes != 0) {
		LOG_ERR("img_DoBrent failed. exiting.");
		goto EXIT;
	}

	if (xMin >= 0) {
		sprintf(strMessage,"Brent found min at gradient scale factor %g. Updating vectors.",xMin);
		TRACE(LOG_MED_PRI,strMessage);
	}
	else {
		iRes =UPDATE_LAMBDA;
		sprintf(strMessage,"img_DoBrent returned minimum for negative scale factor %g.",xMin);
		TRACE(LOG_MED_PRI,strMessage);
		goto EXIT;
	}

	/* scale the img/src by the result */	
	for(i=0; i<g_iSizeSrc; i++) {
		pSrc->pImage[i] += xMin * pdelSrc[i];
	}

	sprintf(strMessage,"Updating source and image.");
	TRACE(LOG_LOW_PRI,strMessage);

	/* condition the source so there are no negative values */
	if(g_Lambda != 0) img_ConditionImage(pSrc,SMALL_VAL);

EXIT:
	TRACE_OUT;
	if (pImgDelta != NULL) {
		lv_free_image_struct(pImgDelta);
		g_pXiImg = NULL;
	}
	return iRes;
}


/***************************
Function:	img_Minimise (based on Num Rec powell's method for conj gradients )
Description:	entry function for conjugate gradient minimisation of a data/lensmodel
				combination. Does a series of 1-d
Arguments:	
Returns:	0: success. non-zero: failure.
****************************/
int	img_Minimise(lv_lensmodel_t *pLens, lv_mapmatrix_t *pMap, lv_image_t *pSource, lv_image_t *pImg,
	lv_image_t *pData, int iMaxIterations, real_t fTolerance, real_t *pFinalChiSq, real_t *pFinalEnt,
	lv_image_t *pPSF, lv_image_t *pNoise, int *iIter) {

	int	j, iRes=0,iHaveBounds=FALSE;
	real_t	gg,gam,fPrevChiSq=0.0,fPrevEnt=0,dgg,fReturnChiSq=0, fRetEntropy=0,fRedChiSqu=0, fPrevLamChi=0,fTargetChiSqu=0;
	real_t	*gSrc=NULL, *hSrc=NULL, *delSrc=NULL, fLower=0, fUpper=FLT_MAX, fJLow=0, fJHigh=0,fMid=0,fJMid=0;
	real_t	d=0,e=0,min1,min2,p,q,r,s,tol1,xm;
	real_t	fSig=0;
	int		iNumImgPix=0,iNumSrcPix=0,iDOF=0;

	TRACE_IN(img_Minimise);

	/* initialise globals, working space and array sizes */
	g_iSizeImg = img_CalcImgSize(pData);
	g_iSizeSrc = img_CalcImgSize(pSource);
	g_pVariance = pNoise;
	g_pData = pData;
	g_pPSF = pPSF;
	g_pMapMatrix = pMap;
	if (g_FixedLambda == -1) {
		g_Lambda = g_Lambda*2.0;
	}
	else {
		g_Lambda = g_FixedLambda;
	}

	g_pTempImg = lv_duplicate_image(pImg);
	g_pTempSrc = lv_duplicate_image(pSource);
	if (g_pTempSrc == NULL || g_pTempImg == NULL) {
		LOG_ERR("failed to duplicate images. No malloc.");
		iRes = 1;
		goto EXIT;
	}

	/* debugging output */
/* create a reverse-projected source map of the inverse pixel variance
	{
		lv_image_t	*pInvVar,*pSrcInvVar;
		int i;

		pInvVar = lv_duplicate_image(pNoise);
		pSrcInvVar = lv_duplicate_image(pSource);
		lv_ZeroRealArray(pSrcInvVar->pImage,img_CalcImgSize(pSrcInvVar));

		for (i=img_CalcImgSize(pInvVar)-1; i>=0; i--) {
			pInvVar->pImage[i] = 1.0/pInvVar->pImage[i];
		}
		img_ReverseProject(pMap, pInvVar, pSrcInvVar->pImage);
		img_DumpImage(pSrcInvVar->pImage,pSrcInvVar->pAxisSize,"revProjVar");

		lv_free_image_struct(pInvVar);
		lv_free_image_struct(pSrcInvVar);
	}
*/

	/* don't forget: non-zero sources and images will be passed as inital values to the
	 * function. This may or may not be what you want */

	/* create the reverse-transposed point spread function and a data image which is
		convolved with the reverse-transposed PSF. Used in calculating the chi-squ gradient */
	if (pPSF != NULL) {
		if (pPSF->pAxisSize[0] != pPSF->pAxisSize[1]) {
			TRACE(LOG_HIGH_PRI,"WARNING: PSF is not a square matrix. Results may be unpredictable");
		}
		g_pRevPSF = lv_duplicate_image(pPSF);
		if (g_pRevPSF == NULL) {
			TRACE(LOG_ERROR_PRI,"Failed to duplicate the PSF. Exiting.");
			iRes =1;
			goto EXIT;
		}
		iRes = img_ReverseTransposeImg(pPSF, g_pRevPSF);
	}

	gSrc = (real_t *)malloc(sizeof(real_t)*g_iSizeSrc);			/* vectors containing previous grad info */
	hSrc = (real_t *)malloc(sizeof(real_t)*g_iSizeSrc);
	delSrc = (real_t *)malloc(sizeof(real_t)*g_iSizeSrc);		/* current gradient vector from chi-squ and ent*/
	if (gSrc == NULL || hSrc == NULL || delSrc == NULL) {
		LOG_ERR("No malloc");
		iRes =1;
		goto EXIT;
	}

	/* copy initial function and gradient values into working space */
	/* this creates an initial "image" if the source is non-zero */
	iRes = img_CalcJ(pImg,pSource,&fPrevChiSq,&fPrevEnt,pNoise,pPSF);
	if (iRes != 0) {
		goto EXIT;
	}
	iRes = img_CalcGradJ(pImg,pSource,delSrc);
	if (iRes != 0) {
		goto EXIT;
	}

	for(j=0; j<g_iSizeSrc; j++) {
		gSrc[j] = -delSrc[j];
		delSrc[j] = hSrc[j] = gSrc[j];
	}

	/* begin conj gradient process. Each iteration does a line
	 *	minimisation of a pseudo-1D line in the direction of
	 *	the gradient pointed to by delSrcFromImg/delSrc
	 */
	for((*iIter)=0; *iIter < iMaxIterations; (*iIter)++) {

		sprintf(strMessage,"Beginning iteration %d. Current Chi-squ: %g, ent: %g, lambda: %g, J: %g",
			(int) *iIter, fPrevChiSq,fPrevEnt,g_Lambda,fPrevChiSq- g_Lambda*fPrevEnt);
		TRACE(LOG_MED_PRI,strMessage);

		/* debugging output */
		if (g_iDebugImgs) {
			img_DumpImage(delSrc,pSource->pAxisSize,"gradSrc");
			img_DumpImage(pSource->pImage,pSource->pAxisSize,"srcIter");
			img_DumpImage(pImg->pImage,pImg->pAxisSize,"imgIter");
		}
		/* do line min. This finds minimum and adjusts source. Doesn't update image */
		iRes = img_LineMin(pImg,pSource,delSrc,&fReturnChiSq,&fRetEntropy);
		if (iRes != 0 && iRes != UPDATE_LAMBDA) {
			goto EXIT;
		}

		/* do we need to increase lambda? */
		/* if LineMin returns an UPDATE_LAMBDA, then the entropy term is dominating
		 	the gradient and making the chi-squ worse. So decrese lamda to let the
			chi-squ have more say in the gradient */
		/* should really only do this if the very first/second iteration can't make any headway */
		if (iRes == UPDATE_LAMBDA) {
				if (*iIter < 2) {
				g_Lambda *= 0.9;
				sprintf(strMessage,"Updating Lambda. Previous: %g, new: %g",g_Lambda/0.9,g_Lambda);
				TRACE(LOG_HIGH_PRI,strMessage);
				/* failsafe... don't let things loop indefinitely for really bad models */
				if (g_Lambda < 1e-4) {
					LOG_ERR("Lambda is imploding! Failing");
					iRes = 1;
					goto EXIT;
				}
				iRes=0;
				goto NEXTSEARCH;
			}
			else {
				iRes=0;
			}
		}

		/* project the new source through the map matrix */
		lv_projectSourceThruMapMatrix(pSource, pImg, g_pMapMatrix);

		/* convolve with the PSF */
		if (g_pPSF != NULL) {
			g_pTempImg = img_ConvolveImgsWithFFT(pImg,g_pPSF,g_pTempImg);
			lv_CopyRealArray(g_pTempImg->pImage,pImg->pImage,g_iSizeImg);
		}

		/* need to update return values etc here. Function could return here
		 * or in a few lines down */
		*pFinalChiSq = fReturnChiSq;
		*pFinalEnt = fRetEntropy;

		/* are we done? */
		if (2.0*fabs(fReturnChiSq-fPrevChiSq) <= fTolerance*(fabs(fReturnChiSq) + fabs(fPrevChiSq) + EPS)) {

			/* not done yet... */
			iNumImgPix = img_CountActiveImgPixels(pImg);
/*
			iNumSrcPix = img_CountNonZeroElements(pSource->pImage,g_iSizeSrc,sqrt(g_imgImgVariance));
*/
			iNumSrcPix = img_CountNonZeroElements(pSource->pImage,g_iSizeSrc,SMALL_VAL*10);
			iDOF = iNumImgPix-iNumSrcPix- pLens->iNumParameters;

			if (iDOF <=0) {
				sprintf(strMessage,"WARNING: iDOF is %d. Setting to 1. img pix: %d, src: %d, params: %d",iDOF,iNumImgPix,iNumSrcPix,pLens->iNumParameters);
				TRACE(LOG_HIGH_PRI,strMessage);
				iDOF=1;
			}

			fSig = sqrt(2.0/iDOF);
			fTargetChiSqu = 1.0 + fSig;
			fRedChiSqu = fReturnChiSq/iDOF;

			/* set this value in case we are doing lambda ajustment with root finding */
			fJMid = fTargetChiSqu-fRedChiSqu;
			sprintf(strMessage,"chi: %g, imgpix: %d, srcpix: %d, params: %d, DOF %d, Target chi-sq: %g, Reduced chi-sq %g",fReturnChiSq,iNumImgPix,iNumSrcPix,pLens->iNumParameters,iDOF,fTargetChiSqu,fRedChiSqu);
			TRACE(LOG_MED_PRI,strMessage);

			/* if the minimisation with current lambda is finished, then check
				to see if the target range of chi-squ has been reached. Adjust
				lambda using brent's method to converge on target
			 */
			if (g_FixedLambda == -1) {
				if (iHaveBounds == FALSE) {

					if( fTargetChiSqu-fRedChiSqu < 0.0) {
						/* this is the upper bound on lambda because the chi-squ isn't good enough. */
						if (g_Lambda < fUpper) {
							fUpper = g_Lambda;
							fJHigh = fTargetChiSqu-fRedChiSqu;
							g_Lambda /= 2;
							sprintf(strMessage,"Found upper bound %g. Target chi-squ minus Reduced: %g",fUpper,fJHigh);
							TRACE(LOG_MED_PRI,strMessage);
						}
					}
					else {
						/* this is the lower bound */
						if (g_Lambda > fLower) {

							/* some models will be so bad they cannot reach the target reduced chi-squ.
								if we hav already increased lambda, but the reduced chi-squ isn't getting
								any better, then just stop */
							if(2.0*fabs(fReturnChiSq-fPrevLamChi) <= fTolerance*(fabs(fReturnChiSq) + fabs(fPrevLamChi) + EPS)){
								sprintf(strMessage,"No improvement after lambda increase. Normal completion after %d iterations. chi-squ: %g, ent: %g, J: %g"
									,(int) *iIter, fReturnChiSq, fRetEntropy, fReturnChiSq- g_Lambda*fRetEntropy);
								TRACE(LOG_HIGH_PRI,strMessage);
								goto EXIT;
							}
							fPrevLamChi=fReturnChiSq;
							fLower = g_Lambda;
							fJLow = fTargetChiSqu-fRedChiSqu;
							g_Lambda *= 2;
							sprintf(strMessage,"Found lower bound %g. Target chi-squ minus Reduced: %g",fLower,fJLow);
							TRACE(LOG_MED_PRI,strMessage);
						}
					}
					/* do we have a bound for the optimal lambda? */
					if (fLower > 0 && fUpper < FLT_MAX) {
						iHaveBounds = TRUE;

						/* fudge the values of the "midpoint" for the first iteration since
							they don't count */
						fMid = fUpper;
						fJMid = fJHigh;
						sprintf(strMessage,"Have bounds %g,%g",fLower,fUpper);
						TRACE(LOG_MED_PRI,strMessage);
					}
				}
				if (iHaveBounds == TRUE) {
					/* have bounds. Use Brent's method to find best lambda. This is a root-finding
						Brent's method, not a minimisation one */
					/* still haven't reached target chi-sq */

					sprintf(strMessage,"fJHigh: %g, fJLow: %g, fJMid: %g, fUpper: %g, fLower: %g, fMid: %g",fJHigh,fJLow,fJMid,fUpper,fLower,fMid);
					TRACE(LOG_MED_PRI,strMessage);
					if((fJMid > 0.0 && fJHigh > 0.0) || (fJMid < 0.0 && fJHigh < 0.0) ) {
						fUpper = fLower;
						fJHigh = fJLow;
						e = d = fJMid - fJLow;
					}
					if (fabs(fJHigh) < fabs(fJMid)) {
						fLower = fMid;
						fMid = fUpper;
						fUpper = fLower;
						fJLow = fJMid;
						fJMid = fJHigh;
						fJHigh = fJLow;
					}
					sprintf(strMessage,"fJHigh: %g, fJLow: %g, fJMid: %g, fUpper: %g, fLower: %g, fMid: %g",fJMid,fJHigh,fJLow,fUpper,fLower,fMid);
					TRACE(LOG_MED_PRI,strMessage);
	
					/* fractional tolerence is the number on the end of next line */
					tol1 = 2.0*EPS*fabs(fMid) + 0.5*0.02;
					xm = 0.5*(fUpper-fMid);

					sprintf(strMessage,"fSig: %g, tol1: %g, xm: %g, fMid: %g",fSig, tol1, xm, fMid);
					TRACE(LOG_MED_PRI,strMessage);

					/* check for doneness */
					if(fabs(xm) <= tol1 || fJMid == 0.0) {
						/* within tolerance. Normal completion */
						sprintf(strMessage,"Normal completion after %d iterations. chi-squ: %g, ent: %g, J: %g, lambda: %g"
							,(int) *iIter, fReturnChiSq, fRetEntropy, fReturnChiSq- g_Lambda*fRetEntropy,g_Lambda);
						TRACE(LOG_HIGH_PRI,strMessage);
						goto EXIT;
					}

					if(fabs(e) >= tol1 && fabs(fJLow) > fabs(fJMid)) {
						s = fJMid/fJLow;
						if (fLower == fUpper) {
							p=2.0*xm*s;
							q=1.0-s;
						}
						else {
							q=fJLow/fJHigh;
							r=fJMid/fJHigh;
							p=s*(2.0 * xm * q * (q-r)-(fMid-fLower)*(r-1.0));
							q=(q-1.0)*(r-1.0)*(s-1.0);
						}
						if (p > 0.0) q = -q;
						p=fabs(p);
						min1 = 3.0*xm*q - fabs(tol1*q);
						min2 = fabs(e*q);
						if(2.0*p < MIN(min1,min2)) {
							e=d;
							d=p/q;
						}
						else {
							d=xm;
							e=d;
						}
					}
					else {
						d=xm;
						e=d;
					}
					fLower=fMid;
					fJLow=fJMid;
					if(fabs(d) > tol1) {
						fMid += d;
					}
					else {
						fMid +=SIGN(tol1,xm);
					}
					g_Lambda = fMid;
					sprintf(strMessage,"Adjusting lambda using Brent. Lower: %g, current: %g, Upper: %g. fJLower: %g, fJMid: %g, fJHigh: %g",fLower,fMid,fUpper,fJLow,fJMid,fJHigh);
					TRACE(LOG_MED_PRI,strMessage);
				}
				/* now reset the source and model image */
/*
*/
				lv_ZeroRealArray(pImg->pImage,g_iSizeImg);
				for(j=0;j<g_iSizeSrc;j++) pSource->pImage[j] = g_DataMean;
/*
				lv_ZeroRealArray(pSource->pImage,g_iSizeSrc);
				img_VecMultByScalar(pImg->pImage,0.5,g_iSizeImg);
				img_VecMultByScalar(pSource->pImage,0.5,g_iSizeSrc);
*/
				iRes = img_CalcGradJ(pImg,pSource,delSrc);
				if (iRes != 0) {
					goto EXIT;
				}
	
				for(j=0; j<g_iSizeSrc; j++) {
					gSrc[j] = -delSrc[j];
					delSrc[j] = hSrc[j] = gSrc[j];
				}
				continue;
			}
			else {
				/* fixed lambda. normal completion */
				sprintf(strMessage,"Normal completion after %d iterations. chi-squ: %g, ent: %g, J: %g, lambda: %g"
					,(int) *iIter, fReturnChiSq, fRetEntropy, fReturnChiSq- g_Lambda*fRetEntropy,g_Lambda);
				TRACE(LOG_HIGH_PRI,strMessage);
				goto EXIT;
			}
		}

NEXTSEARCH:

		/* calculate the new search direction */
		fPrevChiSq = fReturnChiSq;
		fPrevEnt = fRetEntropy;
		lv_ZeroRealArray(delSrc,g_iSizeSrc);
		img_CalcGradJ(pImg,pSource,delSrc);
		if (iRes != 0) {
			goto EXIT;
		}

		dgg = gg = 0.0;

		for (j=0; j < g_iSizeSrc; j++){
			gg += gSrc[j] * gSrc[j];
			dgg += (delSrc[j] + gSrc[j]) * delSrc[j];
		}
		if (gg==0.0) {
			/* unlikely chance that grad is zero then we're done */
			goto EXIT;
		}
		gam = dgg/gg;

		/* update vals for next round */
		for (j=0; j < g_iSizeSrc; j++){
			gSrc[j] = -delSrc[j];
			delSrc[j] = hSrc[j] = gSrc[j] + gam*hSrc[j];
		}
	}

	sprintf(strMessage,"WARNING: Max number of iterations (%d) reached with no convergence",iMaxIterations);
	TRACE(LOG_HIGH_PRI,strMessage);

EXIT:
	if (g_pRevPSF != NULL) {
		lv_free_image_struct(g_pRevPSF);
		g_pRevPSF = NULL;
	}
	if (g_pTempImg != NULL) {
		lv_free_image_struct(g_pTempImg);
		g_pTempImg = NULL;
	}
	if (g_pTempSrc != NULL) {
		lv_free_image_struct(g_pTempSrc);
		g_pTempSrc = NULL;
	}
	if (gSrc!=NULL) free(gSrc);
	if (hSrc !=NULL) free(hSrc);
	if (delSrc != NULL) free(delSrc);
	TRACE_OUT;
	return iRes;
}


/***************************
Image manipulation
****************************/

/***************************
Function:	img_DoConjGradient
Description:	Entry point for this module of code. Commences conjugate gradient minimisation for
				the lensmodel in pLens and (already constructed) mapping matrix pMap.
Arguments:
Returns:	0: success.
****************************/
int	img_DoConjGradient(lv_lensmodel_t *pLens, lv_mapmatrix_t *pMap, lv_image_t *pSource, lv_image_t *pProjectedImage,
	lv_image_t *pMeasuredImage, int iMaxIterations, real_t *pFinalChiSqu, real_t *pFinalEntropy,
	lv_image_t **pBestSource, lv_image_t **pBestImage, lv_image_t *pPSF, lv_image_t *pNoise) {

	real_t	fLastChi = 0.0, fBestChi = FLT_MAX;
	real_t	fNumImgPix = 0.0, fLastEnt=0;
	int		iStatus = 0, iIteration=0;

	TRACE_IN(img_DoConjGradient)

	g_iSizeSrc = img_CalcImgSize(pSource);
	g_iSizeImg = img_CalcImgSize(pProjectedImage);

	fNumImgPix = img_CountActiveImgPixels(pProjectedImage);

	sprintf(strMessage,"starting conj gradient for source (%ld x %ld), image (%ld x %ld). Max iterations %d. Img pixels: %g",
		(long)pSource->pAxisSize[0], (long) pSource->pAxisSize[1], (long) pProjectedImage->pAxisSize[0],
		(long) pProjectedImage->pAxisSize[1], iMaxIterations, fNumImgPix);
	TRACE(LOG_HIGH_PRI,strMessage);

/*
	imgDumpMask(pProjectedImage);
*/

	iStatus = img_Minimise(pLens,pMap,pSource, pProjectedImage, pMeasuredImage, iMaxIterations, MIN_TOL,&fLastChi, &fLastEnt, pPSF, pNoise,&iIteration);
	if (fBestChi > fLastChi && iIteration > 1) {
		fBestChi = fLastChi;
		lv_free_image_struct(*pBestSource);
		lv_free_image_struct(*pBestImage);
		*pBestSource = lv_duplicate_image(pSource);
		*pBestImage = lv_duplicate_image(pProjectedImage);
		sprintf(strMessage,"Best chi-squ: %g. After iteration %d",fLastChi, iIteration);
		TRACE(LOG_MED_PRI,strMessage);
		*pFinalChiSqu = fLastChi;
		*pFinalEntropy = fLastEnt;
	}
	else if (iIteration < 2) {
		/* didn't converge properly. not a valid solution */
		iStatus = 1;
		TRACE(LOG_HIGH_PRI,"img_Minimise returned after < 2 iterations. returning error.");
	}

/* EXIT: */
	TRACE_OUT;
	return iStatus;
}

/***************************
Function:	img_ReverseProject
Description:	projects values in the image plane back into the source plane
				source plane values are a weighted sum of image plane values
Arguments:
Returns:	0:	success, nonzero:	failure
****************************/
int		img_ReverseProject(lv_mapmatrix_t *pMap, lv_image_t *pCurrChiDeriv, real_t *pSourceDelta) {

	mapwt_t		*pWeights;
	lv_axissize_t	i,j;
	size_t		iSizeSrcArr, iSizeImgArr;
	int			iStatus = 0;

	TRACE_IN(img_ReverseProject);

	if (pMap == NULL || pCurrChiDeriv == NULL || pSourceDelta == NULL) {
		LOG_ERR("Input pointer(s) are NULL");
		iStatus =1;
		goto EXIT;
	}

	iSizeImgArr = pMap->dimension[0]*pMap->dimension[1];
	iSizeSrcArr = pMap->dimension[2]*pMap->dimension[3];

	for (j=0; j<iSizeImgArr; j++) {
		pWeights = pMap->array[j];
/*
		if (pWeights != NULL && (pCurrChiDeriv->pMask==NULL || pCurrChiDeriv->pMask[j] == FALSE)) {
*/
		if (pWeights != NULL ) {
			for (i = pMap->min[j] ; i<= pMap->max[j]; i++) {
				if (pMap->pSumSrc[i] > 0) {
                                          pSourceDelta[i] += pWeights[i] * pCurrChiDeriv->pImage[j]/(pMap->pSumSrc[i]);
				}
			}
		}
	}

EXIT:
	TRACE_OUT;
	return iStatus;
}


/***************************
Function:		img_CalcChiSqDeriv
Description:
Arguments:		pImg: real measured data
				pImg2: model image
				pResult: pre-allocated array for chi-squ result.
				iSize: size of pImg and pImg2 arrays.
				pNoise: array of noise data corresponding to the real image data
				pRevPSF: reverse-transposed PSF matrix
Returns:		0: success, nonzero: failure.
****************************/
int	img_CalcChiSqDeriv(lv_image_t *pImg, lv_image_t *pImg2, real_t *pResult, real_t fConstVariance, lv_image_t *pNoise, lv_image_t *pPSF) {

	lv_image_t	*pTempRes=NULL;
	lv_image_t	tempImg;
	int	iStatus = 0;
	size_t	i=0,iSize=0;
	real_t	fVariance = 1.0;

	TRACE_IN(img_CalcChiSqDeriv);

	if (pImg == NULL || pImg2 == NULL || pResult == NULL) {
		iStatus =1;
		goto EXIT;
	}

	iSize = img_CalcImgSize(pImg2);

	/* make working space for temp results */
	pTempRes = lv_duplicate_image(pImg);
	if (pTempRes == NULL) {
		TRACE(LOG_ERROR_PRI,"could not duplicate data image. Exiting.");
		iStatus = 1;
		goto EXIT;
	}

	if (fConstVariance > 0.0) {
		fVariance = fabs(fConstVariance);
	}

	/* fudge to make the convole procedure put the result in pResult */
	tempImg.pImage = pResult;

	for (i = 0; i < iSize; i++) {

		/* image pixels which cannot be affected by the source should not be counted */
/*
*/
		if (pImg2->pMask[i] != FALSE) {
			if (pNoise != NULL) {
				fVariance = pNoise->pImage[i];
			}
			else if (fConstVariance == 0.0) {
				fVariance = fabs(pImg->pImage[i]);
			}

			pTempRes->pImage[i] = 2.0 * (pImg->pImage[i] - pImg2->pImage[i]) / fVariance;
/*
*/
		}
		else {
			pTempRes->pImage[i] = 0.0;
		}
	}

	/* convolve the result with the reverse-transposed PSF */
	if (pPSF != NULL) {
		img_ConvolveImgsWithFFT(pTempRes ,pPSF, &tempImg);
/*
*		img_convolveImages(pTempRes ,pPSF, &tempImg, pMask);
*/
	}
	else {
		lv_CopyRealArray(pTempRes->pImage,pResult,iSize);
	}

EXIT:
	if (pTempRes != NULL) lv_free_image_struct(pTempRes);
	TRACE_OUT;
	return iStatus;
}


/***************************
Function:	img_CalcEntropyDeriv
Description: Calculate the gradient of the entropy of an image
Arguments:	pImg,pResult: image (input) and gradient (output) vectors.
			fArrayTotal: optional image total if it has already been calculated
				to save time. iSize: size of image array
Returns:	0: success. non-zero: failure
****************************/
int	img_CalcEntropyDeriv(lv_image_t *pImg, real_t *pResult) {

	int	iStatus =0,i,iSize;
	double fLogA;

	TRACE_IN(img_CalcEntropyDeriv);

	if (pImg == NULL || pResult == NULL) {
		iStatus =1;
		goto EXIT;
	}
	if (pImg->fDefaultVal <=0) {
		sprintf(strMessage,"Invalid default val: %g",pImg->fDefaultVal);
		LOG_ERR(strMessage);
		iStatus =1;
		goto EXIT;
	}
	fLogA = log(pImg->fDefaultVal);
	iSize = img_CalcImgSize(pImg);

	for (i=0;i<iSize;i++) {
		if (pImg->pImage[i] > 0.0) {
			pResult[i] = -(fLogA - log(pImg->pImage[i]));
		}
		else {
			LOG_ERR("zero/negative source val");
/*
			pResult[i] = -(fLogA - fDefaultVal);
*/
			iStatus=1;
			goto EXIT;
		}
	}

EXIT:
	TRACE_OUT;
	return iStatus;
}


/***************************
Function:	img_MaxArrayVal
Description:	find the max (most positive) value of an array
Arguments:
Returns:	the maximum
****************************/
real_t	img_MaxArrayVal(register real_t *pArr, register size_t iSize) {

	register real_t	fMax = -FLT_MAX;

	if (pArr == NULL) {
		return fMax;
	}

	while (iSize-- >0) {
		if(*pArr > fMax) {
			fMax = *pArr;
		}
		pArr++;
	}

	return fMax;
}


/***************************
Function:	img_CalcArrayTotal
Description:	Calculate the total of an array
Arguments:
Returns:	the total
****************************/
real_t	img_CalcArrayTotal(register real_t *pArr, size_t iSize) {

	register double fTotal =0.0;

	if (pArr == NULL) {
		return 0.0;
	}

	while (iSize-- > 0) {
		fTotal += *(pArr++);
	}
	return (real_t) fTotal;
}

/***************************
Function:	img_CountNonZeroElements
Description:	Count the number of zero elements in an array
				above some threshold
Arguments:
Returns:
****************************/
int	img_CountNonZeroElements(real_t *pArr, register size_t iSize, real_t fThreshold){
	int	iCount=0;
	while (iSize-- > 0) {
		if (pArr[iSize] > fThreshold) {
			iCount++;
		}
	}
	return iCount;
}

/***************************
Function:	img_CalcImgSize
Description: Calc the size of an image array
Arguments:	pImg: pointer to image
Returns: number of pixels in image (not the number of bytes)
*****************************/
int	img_CalcImgSize(lv_image_t *pImg) {
	int	iRes=1,i;

	if (pImg==NULL) return 0;

	for(i=0; i< pImg->iNumAxes; i++) {
		iRes *= pImg->pAxisSize[i];
	}

	return iRes;
}

/***************************
Function:	img_ConditionImage
Description:	Make any negative values in an image into zero.
			Useful for images which will have log operations performed
			on them (such as entropy and derivatives)
Arguments:
Returns:
****************************/
int		img_ConditionImage(lv_image_t *pImg,real_t fDefaultVal) {
	int iSize,iCount=0;
	real_t	*pArr;

	if (pImg == NULL) {
		return 1;
	}

	pArr = pImg->pImage;
	iSize = img_CalcImgSize(pImg);

	while (iSize-- > 0) {
		if (*pArr <= 0) {
			*pArr = fDefaultVal;
			iCount++;
		}
		pArr++;
	}
	return iCount;
}

/***************************
Function:	img_CountActiveImgPixels
Description:
Arguments:
Returns:
****************************/
int	img_CountActiveImgPixels(lv_image_t	*pImg) {

	long	iCount =0, iSize;

	iSize = img_CalcImgSize(pImg);
	while (iSize-- > 0) {
		if (pImg->pMask[iSize] != FALSE) {
			iCount++;
		}
	}
	return iCount;
}

/***************************
Function:	img_CountMultiplyImgdPixels
Description: Count the number of pixels in the image plane which are
			multiply imaged.
Arguments:
Returns:
****************************/
int	img_CountMultiplyImgdPixels(lv_mapmatrix_t *pMap) {
	int	i,iImgSize=0,count=0;

	if (pMap == NULL) return -1;

	iImgSize=pMap->dimension[0]*pMap->dimension[1];

	for (i=0;i<iImgSize;i++) {
		if (pMap->pMultImgPix[i] != FALSE) count++;
	}
	return count;
}


/***************************
Function:	img_FindMultiplyImgdPixels
Description: create an image of the multiply imaged pixels in the
			image plane
Arguments:
Returns:
****************************/
int	img_FindMultiplyImgdPixels(lv_mapmatrix_t *pMap) {
	int	i,j,iSrcSize=0,iImgSize=0,single=0;
	int iStatus=0;

	iSrcSize=pMap->dimension[2]*pMap->dimension[3];
	iImgSize=pMap->dimension[0]*pMap->dimension[1];

	for (i=0;i<iImgSize;i++) {
		single=0;
		pMap->pMultImgPix[i]=FALSE;
		if (pMap->array[i] !=NULL) {
			for (j=pMap->min[i]; j<=pMap->max[i]; j++) {
				if (pMap->array[i][j] > 0 && pMap->pSumSrc[j] > 1.1) {
					single++;
				}
			}
			pMap->pMultImgPix[i]= (single!=0);
		}
	}
/* EXIT:*/
	return iStatus;
}


/***************************
Function:	img_VecScalarProd
Description:	Take the scalar (dot) product of two vectors
Arguments:
Returns:
****************************/
real_t	img_VecScalarProd(register real_t *pArr1, register real_t *pArr2, size_t iSize) {

	register double fTotal = 0.0;

	if (pArr1 == NULL || pArr2 == NULL) {
		return 0.0;
	}

	while (iSize-- > 0) {
		fTotal += *(pArr1++) * *(pArr2++);
	}

	return (real_t) fTotal;
}


/***************************
Function:	img_VecMultByScalar
Description:	Multiply all the values of an array by a scalar.
Arguments:
Returns:
****************************/
int		img_VecMultByScalar(register real_t *pArr, real_t	fScalar, size_t iSize) {

	if (pArr == NULL) {
		return 1;
	}

	while (iSize-- > 0) {
		*(pArr++) *= fScalar;
	}

	return 0;
}


/***************************
Function:
Description:
Arguments:
Returns:
****************************/
int	img_AddArrays(register real_t *pArr1, register real_t *pArr2, register real_t *pResult, size_t iSize) {

	int	iStatus = 0;

	if (pArr1 == NULL || pArr2 == NULL || pResult == NULL) {
		iStatus =1;
		goto EXIT;
	}

	while (iSize-- > 0) {
		*(pResult++) = *(pArr1++) + *(pArr2++);
	}

EXIT:
	return iStatus;
}

/***************************
Function:	img_CalcChiSquared
Description:	Calculate the Chi-Squared difference between two images. Calculated
				as sum( (I(ij) - D(ij)) / sigma(ij) )^2 using calculated variance, or
				as sum( I(ij) - D(ij) )/^2 otherwise
				The variance is calculated as the squre root of the pixel value for Img2
				an error (-ve value) will be returned if a case is encountered where the variance
				is zero and the data difference is non-zero
Arguments:	pImg1- image 1 (calculated image)
			pImg2-	image 2 (real data image)
			pNoise-	image of noise (sigma^2) for the data. If not present, a constant noise value
					(fConstVariance) can be used
Returns:	the value of chi-squre, or -ve number on failure.
****************************/
real_t	img_CalcChiSquared(lv_image_t *pImg1, lv_image_t *pImg2, real_t fConstVariance, lv_image_t *pNoise, bool *pMask) {

	lv_axissize_t	i,iSizeImg=0;
	real_t			fResult =0, fVariance, fDifference, *pImgData1, *pImgData2;

	TRACE_IN(img_CalcChiSquared);

	if (pImg1->pAxisSize[0] != pImg2->pAxisSize[0] || pImg1->pAxisSize[1] != pImg2->pAxisSize[1]) {
		sprintf(strMessage,"Images are not the same size. Can't do Chi-Squared. img1: (%ld x %ld). Img2: (%ld x %ld)",
			(long) pImg1->pAxisSize[0],(long) pImg1->pAxisSize[1],(long) pImg2->pAxisSize[0], (long) pImg2->pAxisSize[1]);
		TRACE(LOG_ERROR_PRI,strMessage);
		fResult = -1.0;
		goto EXIT;
	}

	if (pMask == NULL) {
		sprintf(strMessage,"WARNING: No mask supplied. Results may be wrong.");
		TRACE(LOG_HIGH_PRI,strMessage);
	}

	sprintf(strMessage,"Calculating chi-squared for %ld x %ld images. Const variance (%s): %g",pImg1->pAxisSize[0],
		pImg1->pAxisSize[1],(pNoise == NULL ? "Used" : "Unused"), fConstVariance);
	TRACE(LOG_LOW_PRI,strMessage);

	/* set the variance now if it's always the same */
	fVariance = fConstVariance;
	iSizeImg = img_CalcImgSize(pImg1);

	/* have to do this so pointer arithmetic works. It's faster too. */
	pImgData1 = pImg1->pImage;
	pImgData2 = pImg2->pImage;

	for (i=0; i<iSizeImg; i++) {

		/* only count pixels which are affected by the source */
		if ((pMask == NULL) || pMask[i] != FALSE) {

			/* calc the variance */
			if (pNoise != NULL) {
				fVariance = pNoise->pImage[i];
			}
			else if (fConstVariance == 0.0) {
				fVariance = fabs(pImgData2[i]);
			}
			/* else default variance is set above */

			fDifference = pImgData1[i] - pImgData2[i];

			fResult += (fDifference * fDifference) / fVariance;
		}
	}

EXIT:
	TRACE_OUT;
	return fResult;
}


/***************************
Function:	img_CalcVariance
Description: calculate the variance of an array
Arguments:	pImg: Array to calc variance of
			iArrSize: size of the array
Returns:	variance
****************************/
real_t	img_CalcVariance(real_t *pImg, lv_axissize_t iArrSize) {

	real_t	fTotal = 0.0, fMean = 0.0, fDiff;
	lv_axissize_t	i;

	if (pImg == NULL || iArrSize <=1 ) {
		return 0.0;
	}

	fMean = img_CalcArrayTotal(pImg, iArrSize)/(real_t) iArrSize;
	for (i=0; i< iArrSize; i++) {
		fDiff = *(pImg+i) - fMean;
		fTotal += fDiff * fDiff;
	}

	return fTotal / ( (real_t)iArrSize-1.0);
}


/***************************
Function:	img_CalcImgEntropy
Description:	calculate the entropy of an image. 
Arguments:		pImg: image to calculculate the entropy of
				pResult: result goes in here. 
Returns:	0: Success. nonzero: failure
			puts result in pResult.
****************************/
int		img_CalcImgEntropy(lv_image_t *pImg, real_t *pResult, real_t fSkyVal) {

	real_t	fCurval, fTotalFlux=0.0;
	double	fTotal =0.0;
	int	i,iStatus = 0,iSize;

	TRACE_IN(img_CalcImgEntropy);

	if (pImg == NULL || pResult == NULL) {
		TRACE(LOG_ERROR_PRI,"Null arguments. Exiting.");
		iStatus =1;
		goto EXIT;
	}

	iSize = img_CalcImgSize(pImg);

	/* now calculate the entropy */
	if (fSkyVal ==0) {
		/* calculate the total flux */
		fTotalFlux = img_CalcArrayTotal(pImg->pImage, iSize);
		if (fTotalFlux == 0.0) {
			TRACE(LOG_MED_PRI,"Empty image. No work to do here.");
			goto EXIT;
		}

		for (i=0; i< iSize; i++) {
			if (pImg->pImage[i] > 0.0) {
				fCurval = pImg->pImage[i]/fTotalFlux;
				fTotal += (double) fCurval * log(fCurval);
			}
		}
	}
	else {
		for (i=0; i<iSize;i++) {
			if (pImg->pImage[i] > 0.0) {
				fTotal += pImg->pImage[i]*(log(pImg->pImage[i]/fSkyVal) - 1.0);
			}
		}
	}

EXIT:
	TRACE_OUT;
	*pResult = (real_t) -fTotal;
	return iStatus;
}


/***************************
Function:	img_convolveImages
Description:	convolve an image with a PSF image of some kind. Progressively
			overlaps PSF image with each pixel and multiplies the image pixel
			with the PSF. Generates a new lv_image_t struct as the output if no
			space is provided in pResultImage. Assumes that the PSF is already normalised.
Arguments:	pImg- image to have PSF applied to
			pPSF- the convoltion "kernel" or point-spread-function
			pResultImage- pointer to preallocated space for the result to go in
			pMask-	a "mask" of active pixels in the image. Pixels which have a mask
				value of TRUE should not be considered (i.e. they should be masked out)
Returns:	pointer to convolved image

NOTE:	This function performs a convolution in the strict mathematical sense. It 
		produces the same results as IDL's convol function. From the perspective
		of "overlaying" a point-spread function onto an image, it does not produce
		the expected result unless the psf is symmetric. The resulting image for
		the psf applied to a single pixel is the psf rotated by 180 degrees.

		The action of a psf in astronomy can be thought of as "pushing" the flux from
		a single pixel into surrounding pixels. The convolution function, by definition,
		does not "push" flux. It calculates a pixel's flux by "pulling" flux from
		surrounding pixels into its own. If you draw a picture of this, you can
		see how they are equivalent when the psf is symmetric only.

		The action of "pushing" the flux can be simulated by rotating the PSF by
		180 degrees before being used by this function.
		So users of this function should keep that in mind when coding. If you want
		a true convolution, then just use this funciton as is. If you want to simulate
		the effects of a psf on an image, rotate the psf
		(using img_ReverseTransposeImg) first.

****************************/
lv_image_t *img_convolveImages(lv_image_t *pImg, lv_image_t *pPSF, lv_image_t *pResultImage, bool *pMask) {

	register	lv_axissize_t	k,l;
	register	real_t	fImgVal, fPsfVal, fTotal;
	real_t	*pImgpix, *pPsfpix;
	lv_axissize_t	i,j, iHalfPSFsizex, iHalfPSFsizey, iMinx, iMiny, iMaxx, iMaxy;
	lv_axissize_t	iSizeX, iSizeY, iSizePSFx;
	int				bOddXDim, bOddYDim;

	TRACE_IN(img_convolveImages);

	if (pImg == NULL || pPSF == NULL) {
		TRACE(LOG_ERROR_PRI,"Input pararmeters are NULL. returning.");
		goto EXIT;
	}

	/* calculate constants once */
	iHalfPSFsizex = pPSF->pAxisSize[0]/2;
	iHalfPSFsizey = pPSF->pAxisSize[1]/2;
	iSizePSFx = pPSF->pAxisSize[0];
	bOddXDim = pPSF->pAxisSize[0]%2;
	bOddYDim = pPSF->pAxisSize[1]%2;

	iSizeX = pImg->pAxisSize[0];
	iSizeY = pImg->pAxisSize[1];

	/* duplicate the existing image for results if no space has been provided */
	if (pResultImage == NULL) {
		pResultImage = lv_duplicate_image(pImg);
		if (pResultImage == NULL) {
			TRACE(LOG_ERROR_PRI,"lv_create_image_struct returned NULL pointer. Exiting.");
			goto EXIT;
		}
	}

	/* set pointers to start of image arrays */
	pImgpix = pImg->pImage;
	pPsfpix = pPSF->pImage;

	for (j=0; j < iSizeY ; j++) {

		/* calculate the Y boundaries for PSF overlap */
		iMiny = MAX(0,j-iHalfPSFsizey);
		iMaxy = MIN(iSizeY, j+iHalfPSFsizey + bOddYDim);

		for (i=0; i < iSizeX ; i++) {

			pResultImage->pImage[i + j*iSizeX] = 0.0;

			/* skip pixels which aren't mapped in the unconvolved image */
/*
			if ((pMask!=NULL && pMask[i + j*iSizeX]) || (pImgpix[i + j*iSizeX] == NULL_PIXEL)) {
*/

/*
* THE FOLLOWING LINES ARE BAD FOR THE "PULL" METHOD OF CONVOLUTION! DON'T USE THEM!
*			if (pImgpix[i + j*iSizeX] == NULL_PIXEL) {
*				continue;
*			}
*/

			/* calculate the X boundaries for the PSF overlap */
			iMinx = MAX(0,i-iHalfPSFsizex);
			iMaxx = MIN(iSizeX, i+iHalfPSFsizex + bOddXDim);

			fTotal = 0.0;

			/* overlay the PSF at this image pixel */
			/* l and k loop through the part of the PSF which is not outside the image */
			for (l=iMiny; l < iMaxy; l++) {
				for (k = iMinx; k < iMaxx; k++) {
					fPsfVal = pPsfpix[(k-iMinx) + (l-iMiny) * iSizePSFx];
					fImgVal = pImgpix[k + l * iSizeX];
					if (fImgVal != NULL_PIXEL) {
						fTotal += fImgVal * fPsfVal;
					}
				}
			}
			pResultImage->pImage[i + j*iSizeX] = fTotal;
		}
	}

EXIT:
	TRACE_OUT;
	return pResultImage;
}

/***************************
Function:	img_ConvolveImgsWithFFT
Description:	 convolves an image (pImg) with a kernel (pPSF) using the fast
				fourier transform method. If pResultImage is NULL, it creates
				space for the resulting image, otherwise the result is put in
				pResultImage
Arguments:
Returns:	pointer to lv_image_t struct which contains convolved data
****************************/
lv_image_t	*img_ConvolveImgsWithFFT(lv_image_t *pImg, lv_image_t *pPSF, lv_image_t *pResultImage) {

	static	bool	bDonePlan = FALSE;
	static	rfftwnd_plan	forward_plan, inv_plan;
	static	int	iSizeData[MAX_FT_DIMENSION]={0,0,0,0},iSizeFT[MAX_FT_DIMENSION]={0,0,0,0}, iSizePad[MAX_FT_DIMENSION]={0,0,0,0};
	static	fftw_complex	*pFTImg=NULL, *pFTConv=NULL;
	static	fftw_real	fScale=1;

	int	i,j, iNumElements=1, iStatus=0;
	int	iIsOdd=0, iSizeImg=0;
	fftw_complex	*pTempFFTInfo=NULL;
	fftw_real	*pTemp=NULL;
	lv_image_t	*pTempImg=NULL;

	TRACE_IN(img_ConvolveImgsWithFFT);

	/* unlikely for normal data, but check anyway. */
	if (pImg->iNumAxes > MAX_FT_DIMENSION) {
		sprintf(strMessage,"ERROR: Dimension of image too large: %d. FT won't work.",(int)pImg->iNumAxes);
		LOG_ERR(strMessage);
		return NULL;
	}

	iSizeImg = pImg->pAxisSize[0]*pImg->pAxisSize[1];

	/* create plan and FT of PSF image if it hasn't been done already */
	if (bDonePlan == FALSE) {

		/* set up fftw threads.  Once only.*/
		/*
		if ( fftw_threads_init() !=0)  {
			LOG_ERR("fftw_threads_init failed");
			goto EXIT;
		}
		*/

		/* calculate the sizes of the arrays for the image and FT transform data */
		/* all the axis sizes are the same except the last which is n/2 + 1 */
		/* the data must be padded with zeros to prevent aliasing in the final convolution */
		/* the size of the padding must be the size of the PSF array in both directions */
		for (i=0; i < pImg->iNumAxes; i++) {
			iSizePad[i] = MIN(pImg->pAxisSize[i],pPSF->pAxisSize[i])/2 + 1;
			if (iSizePad[i]%2 == 1) {
				iSizePad[i] += 1;
			}
			iSizeData[i] = MAX(pImg->pAxisSize[i],pPSF->pAxisSize[i]) + iSizePad[i];
			iSizeFT[i] = iSizeData[i];
			fScale /= iSizeData[i];
		}
		iSizeFT[pImg->iNumAxes-1] = iSizeData[pImg->iNumAxes-1]/2 + 1;

		/* first create the execution plan */
		forward_plan = rfftwnd_create_plan(pImg->iNumAxes,iSizeData,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE);
		inv_plan = rfftwnd_create_plan(pImg->iNumAxes,iSizeData,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE);

		sprintf(strMessage,"Created FFT plans for (%dx%dx%dx%d) data",(int) iSizeData[0],(int) iSizeData[1],(int) iSizeData[2],(int) iSizeData[3]);
		TRACE(LOG_MED_PRI,strMessage);

		bDonePlan = TRUE;
	}

	/* calculate the number of elements in the complex FT array */
	for (i=0; i<pImg->iNumAxes;i++) {
		iNumElements *= iSizeFT[i];
	}

	/* calculate the FT of the PSF */
	if (pPSF->ftInfo == NULL) {

		/* pad the PSF image so that its is the same size as the data.
		 * this is necessary for the FFT method. */
		iStatus = img_PadImage(pPSF,iSizeData);
		if (iStatus != 0) {
			TRACE(LOG_ERROR_PRI,"img_PadImage failed. Exiting");
			pResultImage=NULL;
			goto EXIT;
		}
/*
*		if (g_iDebugImgs) img_DumpImage(pPSF->pImage,iSizeData,"padPSF");
*/

		/* now move the image so that it is centred around the zero'th pixel */
		/* for an image this has the effect of chopping the image through the vertical
		and horizontal centre lines, then swapping the NW/SE and SW/NE chunks without
		flipping or rotating the images */

		TRACE(LOG_MED_PRI,"Converting PSF image so that signal is at edges");
		pTemp = malloc(iSizeData[0]*iSizeData[1]*sizeof(fftw_real));
		for (i=0; i<iSizeData[0]*iSizeData[1]; i++) pTemp[i]=0;
		iIsOdd = (iSizeData[0]*iSizeData[1])%2;
		for (j=0; j<iSizeData[1]/2; j++) {
			for (i=0; i< iSizeData[0]/2; i++) {
				pTemp[i+j*iSizeData[0]] = pPSF->pImage[i+iSizeData[0]/2 + (j+iSizeData[1]/2)*iSizeData[0]];
				pTemp[i+iSizeData[0]/2 +iIsOdd + (j+iSizeData[1]/2 + iIsOdd)*iSizeData[0]] = pPSF->pImage[i+j*iSizeData[0]];
				pTemp[i + (j+iSizeData[1]/2+iIsOdd)*iSizeData[0]] = pPSF->pImage[i+iSizeData[0]/2 + j*iSizeData[0]];
				pTemp[i+iSizeData[0]/2+iIsOdd+ j*iSizeData[0]] = pPSF->pImage[i+(j+iSizeData[1]/2)*iSizeData[0]];
			}
		}
/*
*		if (g_iDebugImgs) img_DumpImage(pTemp,iSizeData,"padShiftPSF");
*/

		/* allocate space for the FT data */
		pPSF->ftInfo = (void *) calloc(iNumElements,sizeof(fftw_complex));
		if (pPSF->ftInfo == NULL) {
			TRACE(LOG_ERROR_PRI,"No malloc. Exiting");
			pResultImage=NULL;
			goto EXIT;
		}

		/* create FT of PSF just once */
		TRACE(LOG_MED_PRI,"Calculating FT of PSF");
		/* with threads
		rfftwnd_threads_one_real_to_complex(2,forward_plan,pTemp,pPSF->ftInfo);
		 */
		/* no theads 
		*/
		rfftwnd_one_real_to_complex(forward_plan,pTemp,pPSF->ftInfo);

		free(pTemp);
		pTemp = NULL;
	}

	/* create FT of image data and space for convolution result*/
	if (pFTImg == NULL) {
		pFTImg= (void *) calloc(iNumElements,sizeof(fftw_complex));
	}
	if (pFTConv == NULL) {
		pFTConv= (void *) calloc(iNumElements,sizeof(fftw_complex));
	}
	if (pFTImg == NULL || pFTConv == NULL) {
		TRACE(LOG_ERROR_PRI,"No malloc. Exiting");
		pResultImage=NULL;
		goto EXIT;
	}

	/* duplicate the existing image for results if no space has been provided */
	if (pResultImage == NULL) {
		pResultImage = lv_duplicate_image(pImg);
		if (pResultImage == NULL) {
			TRACE(LOG_ERROR_PRI,"lv_create_image_struct returned NULL pointer. Exiting.");
			goto EXIT;
		}
	}

	/* pad the image then do the FT */
	/* create space for the padded image */
	pTempImg = lv_create_image_struct(pImg->pAxisSize[0],pImg->pAxisSize[1],pImg->iFitsBitPix,pImg->fPixelAngSize);
	if (pTempImg == NULL) {
		TRACE(LOG_ERROR_PRI,"lv_duplicate_image failed. Exiting");
		goto EXIT;
	}

/******* copy image data into temp space and replace "NULL_PIXEL" values with zero
 * attempt workaround for slow FT calculation which has small values in it 
 ******** unnecessary when NULL_PIXEL is 0.0
 *************/
	if (NULL_PIXEL != 0.0) {
		for (i=0;i<iSizeImg; i++) {
			if (pImg->pImage[i] == NULL_PIXEL) {
				pTempImg->pImage[i]=0;
			}
			else {
				pTempImg->pImage[i] = pImg->pImage[i];
			}
		}
	}
	else {
		for (i=0;i<iSizeImg; i++) {
			pTempImg->pImage[i] = pImg->pImage[i];
		}
	}

	iStatus = img_PadImage(pTempImg,iSizeData);
	if (iStatus != 0) {
		TRACE(LOG_ERROR_PRI,"img_PadImage failed for image. Exiting");
		pResultImage = NULL;
		goto EXIT;
	}
	if (g_iDebugImgs) img_DumpImage(pTempImg->pImage,pTempImg->pAxisSize,"imgPad");

	/* calculate the FT of the image */
	/* must copy the image to the data type specified by the FT.
	 * this is inefficient, but can't do much about it if we keep
	 * separate data types for fftw_real and real_t. Do the mult by
	 * fScale here so that we can ignore it a few lines down */

	TRACE(LOG_LOW_PRI,"Calculating FT of image");
	j = pTempImg->pAxisSize[0]*pTempImg->pAxisSize[1];
	pTemp = malloc(sizeof(fftw_real) * j);
	if (pTemp==NULL) {
		LOG_ERR("No malloc to copy image data into fftw_real array");
		pResultImage = NULL;
		goto EXIT;
	}
	for (i=0; i< j; i++) pTemp[i] = pTempImg->pImage[i]* fScale;
	/* with threads
	rfftwnd_threads_one_real_to_complex(2,forward_plan,pTemp,pFTImg);
	*/
	/* no threads
	*/
	rfftwnd_one_real_to_complex(forward_plan,pTemp,pFTImg);

/* debugging. Dump the magnitude of the transformed PSF and image*/
/* Debugging code
*	pTemp = calloc(iNumElements,sizeof(real_t));
*	pTemp2 = calloc(iNumElements,sizeof(real_t));
*	pTempFFTInfo = (fftw_complex *) pPSF->ftInfo;
*	for (i=0; i< iNumElements; i++) {
*		pTemp[i] = sqrt(pTempFFTInfo[i].re*pTempFFTInfo[i].re + pTempFFTInfo[i].im*pTempFFTInfo[i].im);
*	}
*	pTempFFTInfo = (fftw_complex *) pFTImg;
*	for (i=0; i< iNumElements; i++) {
*		pTemp2[i] = sqrt(pTempFFTInfo[i].re*pTempFFTInfo[i].re + pTempFFTInfo[i].im*pTempFFTInfo[i].im);
*	}
*	img_DumpImages(pTemp2,pTemp,(size_t *)iSizeFT,(size_t *)iSizeFT);
*	free(pTemp);
*	free(pTemp2);
*/

	/* now convolve the data in fourier space */
	/* scale has been multiplied already above */
	pTempFFTInfo = (fftw_complex *) pPSF->ftInfo;
	for (i=0; i< iNumElements; i++) {	/* Y values or number of rows */
		pFTConv[i].re = (pFTImg[i].re*pTempFFTInfo[i].re - pFTImg[i].im*pTempFFTInfo[i].im);
		pFTConv[i].im = (pFTImg[i].re*pTempFFTInfo[i].im + pFTImg[i].im*pTempFFTInfo[i].re);
	}

/* debugging code
*	pTemp2 = calloc(iNumElements,sizeof(real_t));
*	for (i=0; i< iNumElements; i++) {
*		pTemp2[i] = sqrt(pFTConv[i].re*pFTConv[i].re + pFTConv[i].im*pFTConv[i].im);
*	}
*	img_DumpImages(pTemp2,NULL,(size_t *)iSizeFT,(size_t *)NULL);
*	free(pTemp2);
*/

	/* invert the FFT to get the convolved images */
	TRACE(LOG_LOW_PRI,"Calculating inverse FT of convolved data");
	/* with threads
	rfftwnd_threads_one_complex_to_real(2,inv_plan,pFTConv,pTemp);
	 */
	/* no threads
	*/
	rfftwnd_one_complex_to_real(inv_plan,pFTConv,pTemp);

	/* this image will be the oversized one. We now need to chop out the
	 * bit we are interested in */
	TRACE(LOG_LOW_PRI,"Extracting real data from padded convolved data");
	for (i=0; i<pImg->pAxisSize[1]; i++) {
		for(j=0; j<pImg->pAxisSize[0]; j++) {
			pResultImage->pImage[j + i*pImg->pAxisSize[0]] = pTemp[(j+iSizePad[0]/2) + iSizeData[0]*(iSizePad[1]/2+i)];
		}
	}

	TRACE(LOG_LOW_PRI,"Freeing temporary image");
EXIT:
	if (pTempImg != NULL) lv_free_image_struct(pTempImg);
	if (pTemp != NULL) free(pTemp);
	TRACE_OUT;
	return pResultImage;
}

/***************************
 * Function:	img_PadImage
 * Description:	Insert an image into a larger one padded with zeros. The original image in
 *				the data structure is free'd.
 * Arguments:	pImg: Image to insert into bigger one
 * 				pDimension: array of dimension sizes to be expanded to
 * Returns:		0: success, nonzero: failure
 *****************************/
int img_PadImage(lv_image_t *pImg, int *pDimension) {
	int iStatus = 0;
	lv_axissize_t   iSize=0, i,j, iXDiff, iYDiff;
	real_t  *pNewImg = NULL;
	bool    bIsSameSize = TRUE;

	TRACE_IN(img_PadImage);
	if (pImg->iNumAxes != 2) {
		TRACE(LOG_ERROR_PRI,"Only works for 2-d images at the moment.");
		iStatus =1;
		goto EXIT;
	}

	for (i=0; i<pImg->iNumAxes; i++) {
		if (pDimension[i] != pImg->pAxisSize[i]) {
			bIsSameSize = FALSE;
			break;
		}
	}

	if(bIsSameSize) {
		TRACE(LOG_LOW_PRI,"Sizes are same. Doing nothing.");
		goto EXIT;
	}

	sprintf(strMessage,"Padding image. Previous size: %dx%d. New size %dx%d.",(int)pImg->pAxisSize[0],
					(int)pImg->pAxisSize[1],(int)pDimension[0],(int)pDimension[1]);
	TRACE(LOG_LOW_PRI,strMessage);

	iSize = pDimension[0] * pDimension[1];
	pNewImg = malloc(sizeof(real_t) * iSize);
	if (pNewImg == NULL) {
		TRACE(LOG_ERROR_PRI,"No malloc");
		iStatus=1;
		goto EXIT;
	}

	lv_ZeroRealArray(pNewImg,iSize);
	iXDiff = (pDimension[0] - pImg->pAxisSize[0])/2 + (pDimension[0] - pImg->pAxisSize[0])%2;
	iYDiff = (pDimension[1] - pImg->pAxisSize[1])/2 + (pDimension[0] - pImg->pAxisSize[0])%2;
	for (j=0; j<pImg->pAxisSize[1]; j++) {
		for (i=0; i< pImg->pAxisSize[0]; i++) {
			pNewImg[i + pDimension[0]*(j+iYDiff) + iXDiff] = pImg->pImage[i + j*pImg->pAxisSize[0]];
		}
	}
	free(pImg->pImage);
	pImg->pImage = pNewImg;
	for (i=0; i<pImg->iNumAxes; i++) {
		pImg->pAxisSize[i] = pDimension[i];
	}

EXIT:
	TRACE_OUT;
	return iStatus;
}


/***************************
Function:	img_ReverseTransposeImg
Description:	Flip an image (matrix, whatever) along both diagonals.
				Equivalent to a rotation by 180 degrees
Arguments:		pImg: Image to be reverse-transposed
				pResult:	destination. Must be already allocated
Returns:		0- success. Non-zero: failure
****************************/
int		img_ReverseTransposeImg(lv_image_t *pImg, lv_image_t *pResult) {
	size_t	i,iSize;
	int		iResult=0;

	TRACE_IN(img_ReverseTransposeImg);
	if (pImg == NULL || pResult == NULL) {
		sprintf(strMessage,"Arguments are NULL. Exiting.");
		TRACE(LOG_ERROR_PRI,strMessage);
		iResult=1;
		goto EXIT;
	}

	iSize = img_CalcImgSize(pImg);

	for (i=0; i<iSize; i++) {
		pResult->pImage[iSize-i-1] = pImg->pImage[i];
	}

EXIT:
	TRACE_OUT;
	return iResult;
}


/***************************
Function:
Description:
Arguments:
Returns:
****************************/
lv_image_t	*img_CalcImgPlaneInvMag(lv_mapmatrix_t *pMap, real_t fSrcPixAngSize, real_t fImgPixAngSize) {
	int	i,iSizeImg=0;
	lv_image_t  *pMag = NULL;
	real_t	fRatio=1.;

	TRACE_IN(img_CalcImgPlaneInvMag);

	pMag = lv_create_image_struct(pMap->dimension[0],pMap->dimension[1],IMGTYPE,fImgPixAngSize);
	if (pMag == NULL) {
		TRACE(LOG_ERROR_PRI,"lv_create_image_struct failed");
		goto EXIT;
	}

	iSizeImg = pMap->dimension[0]*pMap->dimension[1];
	fRatio = (fSrcPixAngSize*fSrcPixAngSize)/(fImgPixAngSize*fImgPixAngSize);
	for(i=0; i< iSizeImg; i++) {
		if (pMap->pSumImg[i] > 0) {
			pMag->pImage[i]=(pMap->pSumImg[i] * fRatio);
		}
		else {
			pMag->pImage[i]=1.;
		}
	}

EXIT:
	TRACE_OUT;
	return pMag;
}

/***************************
Function:	img_CalcSrcPlaneMag
Description:	Calculate the magnification for each pixel in the source
				plane. Can be used to generate caustics etc.
Arguments:	pMap:	magnification mapping matrix
Returns:	pointer to image struct, or NULL on failure
****************************/
lv_image_t	*img_CalcSrcPlaneMag(lv_mapmatrix_t *pMap, real_t fSrcPixAngSize, real_t fImgPixAngSize) {

	lv_image_t	*pSrcMag = NULL;
	real_t	fPixGridSizeRatio = 0.0, fSrcTotal;
	mapwt_t	*pWeights;
	lv_axissize_t   i,j,iSizeSrc=0, iSizeImg=0;

	TRACE_IN(img_CalcSrcPlaneMag);

	/* create space for the resulting image */
	pSrcMag = lv_create_image_struct(pMap->dimension[2],pMap->dimension[3],IMGTYPE,fSrcPixAngSize);
	if (pSrcMag == NULL) {
		TRACE(LOG_ERROR_PRI,"lv_create_image_struct returned NULL.");
		goto EXIT;
	}

	fPixGridSizeRatio = (fImgPixAngSize*fImgPixAngSize)/(fSrcPixAngSize*fSrcPixAngSize);
	iSizeSrc = img_CalcImgSize(pSrcMag);
	iSizeImg = pMap->dimension[0]*pMap->dimension[1];

	lv_ZeroRealArray(pSrcMag->pImage, iSizeSrc);


	/* for each source pixel,
	 * find all image pixels which map to it & sum their area
	 * find total area of in source plane of the mapped image pixels
	 * inv mag = (src pix size/src plane total)/image area total
	 */

	for (i=0; i<iSizeSrc; i++) {
		fSrcTotal = 0.0;
		for(j=0; j<iSizeImg; j++) {
			pWeights = pMap->array[j];
			if (pWeights == NULL) {
				continue;
			}
			if (pWeights[i] > 0.0) {
				fSrcTotal += pWeights[i]/pMap->pSumImg[j];
			}
		}
		pSrcMag->pImage[i] = fPixGridSizeRatio*(fSrcTotal);
	}

EXIT:
	TRACE_OUT;
	return pSrcMag;
}


/**************************
debugging code
* useful if you're stepping through the code in a debugger and want to see
* what the images at intermediate steps look like.
***************************/
void	img_DumpImage(real_t *pData, lv_axissize_t *pSizeArr,char *pName) {
	static int	count=0;
	char	strFileName[80];
	lv_image_t  outImg;
	lv_axissize_t	Axis[2];

	if (pName ==NULL) {
		sprintf(strFileName,"imgdump%03d.fits",count);
	}
	else {
		sprintf(strFileName,"%s%03d.fits",pName,count);
	}
	outImg.iNumAxes = 2;
	outImg.iFitsBitPix =IMGTYPE;
	outImg.fPixelAngSize = g_PixelResn;
	outImg.iTableDataType =  lv_decode_FITS_table_data_size(IMGTYPE);
	outImg.iPixelDepth = lv_decode_FITS_pix_byte_size(IMGTYPE);
	outImg.pAxisSize = Axis;
	outImg.pImage = pData;

	if (pData != NULL) {
		count++;
		Axis[0] = pSizeArr[1];
		Axis[1] = pSizeArr[0];
		lv_write_image_to_file(&outImg, strFileName, TRUE);
	}
}

void	imgDumpMask(lv_image_t *pImg) {
	lv_image_t	*pTemp=NULL;
	float	*pMask=NULL;
	int i,iSize=0;

	pTemp = lv_create_image_struct(pImg->pAxisSize[0],pImg->pAxisSize[1],FLOAT_IMG,pImg->fPixelAngSize);
	if (pTemp ==NULL) return;
	iSize = img_CalcImgSize(pImg);
	/* force the pointer to be pointer to char, since we have defined this image
	 * as BYTE_IMG */
	pMask = (float *) pTemp->pImage;
	for (i=0; i< iSize; i++) {
		pMask[i] = (float) pImg->pMask[i];
	}
	lv_write_image_to_file(pTemp, "imgmask.fits", TRUE);
	lv_free_image_struct(pTemp);
}

