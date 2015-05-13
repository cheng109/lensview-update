/******************************
lv_mem.c:     implementation of the Skilling & Bryan (1984) MEM for
			lensed image deconvolution
Copyright (C) 2006. Randall Wayth.
$Id: lv_mem.c,v 1.10 2008/10/31 20:43:10 rwayth Exp rwayth $
*******************************/
#include	<unistd.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>
#include	"fitsio.h"
#include	"log.h"
#include	"common.h"
#include	"lv_common.h"
#include	"lv_image.h"
#include	"gsl/gsl_cblas.h"
#include	"gsl/gsl_eigen.h"
#include	"gsl/gsl_matrix.h"
#include	"lv_mem.h"

#define	NUM_SEARCH_DIR	3
#define	GOLD_RATIO	1.618034
#define	ALPH_SRCH_MAX	1000
#define	ALPH_SRCH_MIN	1e-3
#define	PENALTY_MAX	1e6
#define	SMALL_EVAL	5e-5   /* approx floating point rounding error */

/****************************
 * public global variables 
 * ***************************/
int	g_bUseMultImgPixOnly = FALSE;

/****************************
private functions
*****************************/
static	int mem_CalcGradGradC(real_t *pGradGradC, lv_image_t *pNoise, lv_mapmatrix_t *pMap, lv_image_t *pPSF, lv_image_t *pRevPSF);
static	int	mem_CalcSearchDirs(lv_mapmatrix_t *pMap, lv_image_t *pData, lv_image_t *pImg, lv_image_t *pSrc, lv_image_t *pNoise, real_t *pSearchVector[NUM_SEARCH_DIR], real_t *pGradGradC, lv_image_t *pPSF,real_t *pGradS, real_t *pGradC, real_t *pTEST);
static int	mem_CalcLinIndepModels(real_t *pSearchVectors[], real_t *pLinIndVectors[], real_t *pGradGradC, real_t pM_Eigenvales[], int pUseVec[NUM_SEARCH_DIR]);
static int	mem_CalcQ(lv_mapmatrix_t *pMap,lv_image_t *pSrc, lv_image_t *pImg, real_t *pResultChiSqu, real_t *pResEntropy, lv_image_t *pPSF);

/****************************
private global variables
*****************************/
static  int     g_iSizeImg=0;
static  int     g_iSizeSrc=0;
static	lv_image_t	*g_pSrc=NULL;
static	lv_image_t	*g_pTempImg=NULL;
static	lv_image_t	*g_pData=NULL;
static	lv_image_t	*g_pRevPSF=NULL;
static	lv_image_t	*g_pVariance=NULL;

/***************************
Function:   mem_DoMEM
Description:    Entry point for this module of code.
Arguments:
Returns:    0: success.
****************************/
int mem_DoMEM(lv_lensmodel_t *pLens, lv_mapmatrix_t *pMap, lv_image_t *pSource, lv_image_t *pProjectedImage,
		lv_image_t *pMeasuredImage, int iMaxIterations, real_t *pFinalChiSqu, real_t *pFinalEntropy,
		lv_image_t **pResultSource, lv_image_t **pResultImage, lv_image_t *pPSFin, lv_image_t *pNoise, real_t fSrcDefaultVal,
		real_t  *pFinalQ) {
	real_t  fChi = 0.0, fLastChi = 0,fEnt=0, fTEST=0, fQ=0;
	real_t	fX_mu[NUM_SEARCH_DIR], fBestX_mu[NUM_SEARCH_DIR], *pLinIndVector[NUM_SEARCH_DIR];
	real_t	*pSearchVector[NUM_SEARCH_DIR], fS_mu[NUM_SEARCH_DIR], fC_mu[NUM_SEARCH_DIR];
	real_t	*pGradGradC = NULL, *pGradS=NULL, *pGradC=NULL, fM_Eigenvals[NUM_SEARCH_DIR];
	int		iUseSearchVec[NUM_SEARCH_DIR]={1,1,1},iNoImprov=0, count=0,loopcount=0;
    int     iStatus = 0, iIteration=0,i,j, iNumSrcPix=0,iNumImgPix=0,iNumMultImgPix=0;
	double	fPenalty=0, fAlpha=0, fAlphaMin=0, fCaim=0, fLo_squ=0, fC_0=0,fC_p=0,fC=0;
	double	fCIntermediateAim=0,fCmin=0,fL_squ=0, fLastAlpha=0,fAlphaMax=ALPH_SRCH_MAX,fAlphaSearchMin=0;
	double	fPMax = PENALTY_MAX, fPMin = 0;
	lv_image_t *pPSF=NULL;
	bool	bAlphaChopDone=FALSE, bPChopDone=FALSE, bAlphaChopSuccess=FALSE;

    TRACE_IN(mem_DoMEM)

	/* initialisation. Set values which do not change throughout the process
		and useful pointers to commonly used stuff */
    g_iSizeSrc = img_CalcImgSize(pSource);
    g_iSizeImg = img_CalcImgSize(pProjectedImage);
	pSource->fDefaultVal = fSrcDefaultVal;
	g_pSrc = lv_duplicate_image(pSource); 		/* current working source */
	g_pData = pMeasuredImage;					/* the data */
	g_pVariance = pNoise;
	/* the PSFs are modified by the FFT process, so make copies to so they aren't damaged */
	pPSF = lv_duplicate_image(pPSFin);
	g_pRevPSF = lv_duplicate_image(pPSF);
	/* note: must use pProjectedImage so that the mask is duplicated */
	g_pTempImg = lv_duplicate_image(pProjectedImage);
	pGradGradC = calloc(g_iSizeSrc,sizeof(real_t));
	for (i=0; i<NUM_SEARCH_DIR; i++) {
		pSearchVector[i]=NULL;
		pLinIndVector[i]=NULL;
	}
	if (g_pSrc == NULL || (pPSFin != NULL && g_pRevPSF == NULL) || g_pTempImg == NULL || pGradGradC==NULL) {
		LOG_ERR("Failed to duplicate PSF or alloc temp img working space.");
		iStatus=1;
		goto EXIT;
	}
	pGradS = malloc(g_iSizeSrc * sizeof(real_t));
	pGradC = malloc(g_iSizeSrc * sizeof(real_t));
	if (pGradS == NULL || pGradC==NULL) {
		LOG_ERR("No malloc for pGradS or pGradC");
		iStatus=-1;
		goto EXIT;
	}

	img_ReverseTransposeImg(pPSF, g_pRevPSF);
	/*
	*	if (g_iDebugImgs) img_DumpImage(g_pRevPSF->pImage,g_pRevPSF->pAxisSize,"revPSF");
	*/

    iNumImgPix = img_CountActiveImgPixels(pProjectedImage);
	iStatus = img_FindMultiplyImgdPixels(pMap);
	if (iStatus != 0) goto EXIT;
    iNumMultImgPix = img_CountMultiplyImgdPixels(pMap);
	for (i=0; i<g_iSizeImg; i++ ) {
		if (pMap->pMultImgPix[i]==TRUE && pProjectedImage->pMask[i]==TRUE) count++;
	}
	sprintf(strMessage,"total non-masked image pixels: %d, number of multiply imaged pixels: %d, number overlap %d. Using multiply imaged pixels only: %d",iNumImgPix,iNumMultImgPix,count,g_bUseMultImgPixOnly);
	TRACE(LOG_HIGH_PRI,strMessage);

	/* make space for the search direction images */
	for(i=0; i<NUM_SEARCH_DIR; i++) {
		pSearchVector[i] = malloc(sizeof(real_t)*g_iSizeSrc);
		pLinIndVector[i] = malloc(sizeof(real_t)*g_iSizeSrc);
		if (pSearchVector[i] == NULL || pLinIndVector[i]==NULL) {
			LOG_ERR("Cannot malloc space for search direction images");
			iStatus=1;
			goto EXIT;
		}
	}

	/* initialise working results to zero */
	for (i=0; i<NUM_SEARCH_DIR; i++) {
		fX_mu[i]=0;
		fBestX_mu[i]=0;
		fS_mu[i]=0;
		fC_mu[i]=0;
	}

    sprintf(strMessage,"starting MEM. Max iterations %d. Src default value: %g",iMaxIterations,fSrcDefaultVal);
    TRACE(LOG_HIGH_PRI,strMessage);

	/* set up the initial source with default values.*/
	for (i=0;i<g_iSizeSrc; i++) {
		if(g_pSrc->pImage[i]<=0) g_pSrc->pImage[i] = g_pSrc->fDefaultVal*0.99;
	}

	/* calculate grad(gradC) which is constant throughout the process */
	iStatus = mem_CalcGradGradC(pGradGradC, pNoise, pMap, pPSF,g_pRevPSF);
    /* experiment...
     *	for(i=0; i< g_iSizeSrc; i++) pGradGradC[i]=0;
     */
    
	if (g_iDebugImgs) img_DumpImage(pGradGradC,g_pSrc->pAxisSize,"gradGradCsrc");


	/* create initial projection of source and calculate initial Ent & Chisqu */
	iStatus = mem_CalcQ(pMap,g_pSrc,pProjectedImage,&fChi,&fEnt,pPSF);
	if (g_iDebugImgs) img_DumpImage(g_pTempImg->pImage,pProjectedImage->pAxisSize,"initProj_noPSF");
	if (g_iDebugImgs) img_DumpImage(pProjectedImage->pImage,pProjectedImage->pAxisSize,"initProj_PSF");
	fQ = fEnt*fAlpha - fChi;

	if (iStatus != 0) {
		LOG_ERR("Problems calcing initial Q");
		goto EXIT;
	}

	/* start iterating. Keep going until we're not improving the result or we reach
		the max iterations or the target chisqu is reached */
	for (iIteration=1;iIteration <=iMaxIterations; iIteration++) {

		sprintf(strMessage,"Beginning iteration %d",iIteration);
		TRACE(LOG_MED_PRI,strMessage);

		/* calculate the current DOF used in the source */
		iNumSrcPix = img_CountNonZeroElements(g_pSrc->pImage,g_iSizeSrc,g_pSrc->fDefaultVal*4);

		/* calculate the current search directions */
		/* assumes that pProjectedImage is a current projection of the
			working source */
		iStatus = mem_CalcSearchDirs(pMap, pMeasuredImage, pProjectedImage, g_pSrc, pNoise, pSearchVector,pGradGradC, g_pRevPSF,pGradS,pGradC,&fTEST);
		if(iStatus != 0) {
			goto EXIT;
		}

	/* debugging output
	*/
		if (g_iDebugImgs) {
			img_DumpImage(g_pSrc->pImage,pSource->pAxisSize,"srcIter");
			img_DumpImage(pProjectedImage->pImage,pProjectedImage->pAxisSize,"imgIterS");
		 	img_DumpImage(pGradS,pSource->pAxisSize,"gradS");
			img_DumpImage(pGradC,pSource->pAxisSize,"gradC");
			img_DumpImage(pSearchVector[0],pSource->pAxisSize,"e1");
			img_DumpImage(pSearchVector[1],pSource->pAxisSize,"e2");
			img_DumpImage(pSearchVector[2],pSource->pAxisSize,"e3");
		}

		/* calculate the diagonalised g, linearly independent search directions
			and eigenvalues of M */
		iStatus = mem_CalcLinIndepModels(pSearchVector, pLinIndVector, pGradGradC, fM_Eigenvals, iUseSearchVec);
		if (iStatus != 0) {
			sprintf(strMessage,"Problems calculating diagonalised models. Result: %d",iStatus);
			LOG_ERR(strMessage);
			goto EXIT;
		}
        sprintf(strMessage,"Using search vectors: (%d,%d,%d)\n",iUseSearchVec[0],iUseSearchVec[1],iUseSearchVec[2]);
        TRACE(LOG_MED_PRI,strMessage);

		if (g_iDebugImgs) {
			img_DumpImage(pLinIndVector[0],pSource->pAxisSize,"e1_new");
			img_DumpImage(pLinIndVector[1],pSource->pAxisSize,"e2_new");
			img_DumpImage(pLinIndVector[2],pSource->pAxisSize,"e3_new");
		}

		/* set initial alpha, Chi targets, dist targets etc. */
		fC_0 = fChi;
		fAlphaMin=0;
		fCaim = iNumImgPix - iNumSrcPix - pLens->iNumParameters;
		if (fCaim <0 ) fCaim=0;
		fCaim += sqrt(2.0*fCaim);
		if (g_fTargetChiSqu >= 0) fCaim=g_fTargetChiSqu;
		fCmin = fC_0;
		for (i=0;  i<NUM_SEARCH_DIR; i++) {
            if (iUseSearchVec[i]) {
                fS_mu[i] = img_VecScalarProd(pLinIndVector[i],pGradS,g_iSizeSrc);
                fC_mu[i] = img_VecScalarProd(pLinIndVector[i],pGradC,g_iSizeSrc);
                fCmin -= 0.5*fC_mu[i]*fC_mu[i]/fM_Eigenvals[i];
            }
            else {
                fS_mu[i] = 0;
                fC_mu[i] = 0;
            }
			sprintf(strMessage,"S_%d: %g, C_%d: %g",i,fS_mu[i], i,fC_mu[i]);
			TRACE(LOG_MED_PRI,strMessage);
		}
		for (i=0; i<NUM_SEARCH_DIR; i++) fBestX_mu[i] = 0;
		fCIntermediateAim = MAX(0.67*fCmin+0.33*fC_0,fCaim);
		fPenalty=0;
		fPMin =0;
		fPMax = PENALTY_MAX;
		bPChopDone = FALSE;
		/* magic number multiplier in next line should be between 0.1 and 0.5 */
		fLo_squ = 0.2 * img_CalcArrayTotal(g_pSrc->pImage,g_iSizeSrc);

		sprintf(strMessage,"MEM loop %03d. C_o: %g, Caim: %g, Cmin: %g, Caim (intermediate): %g. Ent: %g. Lo_Squ: %g, TEST: %g, img pix: %d, src pix: %d",iIteration,fC_0,fCaim,fCmin, fCIntermediateAim,fEnt,fLo_squ,fTEST,iNumImgPix,iNumSrcPix);
		TRACE(LOG_MED_PRI,strMessage);

		if (fC_0 <= fCIntermediateAim) break; 

		/* start main loop. */
		loopcount=0;
		do {	/* while P not zero or P chop not done */

			/* initialise alpha. Need current chi-squ, entropy */
			bAlphaChopDone = FALSE;
			bAlphaChopSuccess = FALSE;
			fAlphaSearchMin = ALPH_SRCH_MIN;
			fAlphaMax = ALPH_SRCH_MAX;
			fLastAlpha=0;
			bPChopDone = FALSE;
			if (fEnt > 0) {
				fAlpha = (fChi-fCIntermediateAim)/fEnt;
			}
			else {
				fAlpha=0;
			}
			if (fAlpha < 0) {
				fAlpha=0;
			}

			sprintf(strMessage,"Beginning alpha-chop with chi: %g, ent: %g, alpha: %g, P: %g",fChi,fEnt,fAlpha,fPenalty);
			TRACE(LOG_MED_PRI,strMessage);

			do {	 /* while alpha chop not done */

				/* calculate x_mu, C, C_p, l^2 */
				loopcount++;
				fC_p = fC_0;
				fC = fC_0;
				fL_squ =0;
				for (i=0; i<NUM_SEARCH_DIR; i++) {
                    fX_mu[i]=0;
                    if(iUseSearchVec[i]) {
					    fX_mu[i] = (fAlpha*fS_mu[i] - fC_mu[i])/(fPenalty + fAlpha + fM_Eigenvals[i]);
                        fC_p += fX_mu[i]*fC_mu[i] + 0.5*(fPenalty + fM_Eigenvals[i])*fX_mu[i]*fX_mu[i];
                        fC += fX_mu[i]*fC_mu[i] + 0.5*fM_Eigenvals[i]*fX_mu[i]*fX_mu[i];
                        fL_squ += fX_mu[i]*fX_mu[i];
                    }
				}

				sprintf(strMessage,"Alpha-chop iteration. alpha,min,max: %g,%g,%g last alpha: %g, C: %g, C_p: %g, LSqu: %g, x0: %g, x1: %g, x2: %g",fAlpha,fAlphaSearchMin,fAlphaMax,fLastAlpha,fC,fC_p,fL_squ,fX_mu[0],fX_mu[1],fX_mu[2]);
				TRACE(LOG_LOW_PRI,strMessage);

				/* chop alpha. Expand geometrically if the current movement is in the same
					direction as the previous, otherwise do a sort of bisection */
				if (fC_p > fC_0 ) {
					/* alpha too large */
					fAlphaMax = fAlpha;
					fLastAlpha = fAlpha;
					TRACE(LOG_LOW_PRI,"alpha too large...");
				}
				else if(fC < fCIntermediateAim || fL_squ >= fLo_squ){
					/* alpha too small */
					fAlphaSearchMin = fAlpha;
					fLastAlpha = fAlpha;
					TRACE(LOG_LOW_PRI,"alpha too small...");
				}
				else {
					/* acceptable result. Save and adjust alpha.
					 * we still want to get fC as close to fCIntermediateAim
					 * as we can */
					for (i=0; i < NUM_SEARCH_DIR; i++) {
						fBestX_mu[i] = fX_mu[i];
					} 
					sprintf(strMessage,"Saving Xmus for alpha %g... fC: %g. XMus: %g, %g, %g",fLastAlpha,fC,fX_mu[0],fX_mu[1],fX_mu[2]);
					TRACE(LOG_MED_PRI,strMessage);
					if (fC < fCIntermediateAim) {
						fAlphaSearchMin = fAlpha;
					} else {
						fAlphaMax = fAlpha;
					}
					fAlphaMax = fAlpha;
					fLastAlpha = fAlpha;
				}
				fAlpha = (fAlphaMax + fAlphaSearchMin)/2;

				/* have we reached our intermediate aim? */
/*
				if(fCIntermediateAim <= fC && fC <= fC_0 && fC_p < fC_0 || fabs((fAlphaMax-fAlphaSearchMin)/fAlphaMax) < 0.05 ) 
*/
				/* keep the alpha-chop going until the uncertainty in alpha is small. */
				/* when alpha-chop terminates, we should automatically have reached 
				 the intermediate aim */
				/*
				if( fabs((fAlphaMax-fAlphaSearchMin)/fAlphaSearchMin) < 1e-5 || fabs(fAlphaMax-fAlphaSearchMin)< 1e-4) {
				*/
				if( fabs(fAlphaMax-fAlphaSearchMin) < ALPH_SRCH_MIN) {
					bAlphaChopDone=TRUE;
				}

			} while ( !bAlphaChopDone);

			/* restore successful alpha */
			fAlpha = fLastAlpha;

			sprintf(strMessage,"Alpha-chop done. alpha: %g, C: %g, C_p: %g, Lsqu: %g",fAlpha,fC,fC_p,fL_squ);
			TRACE(LOG_MED_PRI,strMessage);

			/* a successful alpha-chop obtains a result where
				Caim(intermediate) <= C <= Cp <= C_o */
			if (fL_squ < fLo_squ && (fCIntermediateAim <= fC) && (fC_p <= fC_0)) {
				bAlphaChopSuccess=TRUE;

				/* if alpha-chop successful but P-chop isn't done, decrese P */
				if(fPenalty>0 && !bPChopDone) {
					/* decrease P */
					fPMax = fPenalty;
					fPenalty = (fPMax+fPMin)/2;
				}
			}
			else {
				/* alpha-chop not good, so increase P */
				if (fPenalty ==0) {
					/* 1st increase. set to 1/Lo^2 on dimensional grounds */
					fPenalty = 1.0/fLo_squ;
				}
				else {
					fPMin = fPenalty;
					if (fPMax == PENALTY_MAX) {
						/* geometric increase */
						fPenalty *= GOLD_RATIO;
					}
					else {
						fPenalty = (fPMax+fPMin)/2;
					}
				}
			}
			if (fPenalty > 0 && fabs((fPMax-fPMin)/fPenalty) < 0.01) {
				bPChopDone=TRUE;
			}

			sprintf(strMessage,"P: %g (min/max): %g/%g, P chop done : %d. alphachop success: %d",fPenalty,fPMin,fPMax,bPChopDone,bAlphaChopSuccess);
			TRACE(LOG_MED_PRI,strMessage);

			/* for some unknown reasons, the chops get stuck.
			 * Abort if we are spinning in an infinite loop */
			if (loopcount >=10000) {
				LOG_ERR("loopcount large. Aborting chops");
			}

		}	while ( loopcount < 10000 && (!bAlphaChopSuccess || (fPenalty>1.0e-6 && !bPChopDone)) );
			
		/* update best results if we have one */
		sprintf(strMessage,"Updating source with best X_mus: %g, %g, %g",fBestX_mu[0],fBestX_mu[1],fBestX_mu[2]);
		TRACE(LOG_MED_PRI,strMessage);

		for (j=0; j< g_iSizeSrc; j++) {
			for(i=0;i<NUM_SEARCH_DIR; i++) {
				g_pSrc->pImage[j] += fBestX_mu[i]*pLinIndVector[i][j];
			}
		}

		/* condition the source so there are no negative values */
		i = img_ConditionImage(g_pSrc,g_pSrc->fDefaultVal);
		sprintf(strMessage,"Fixed %d zero/negative source vals",i);
		TRACE(LOG_MED_PRI,strMessage);

		/* re-project source and calc Entropy and Chisqu */
		iStatus = mem_CalcQ(pMap,g_pSrc,pProjectedImage,&fChi,&fEnt,pPSF);
		fQ = fEnt*fAlpha - fChi;
		if (iStatus != 0) {
			LOG_ERR("Problems calcing intermediate Q");
			goto EXIT;
		}

		sprintf(strMessage,"done iteration. chisqu: %g, ent: %g",fChi,fEnt);
		TRACE(LOG_MED_PRI,strMessage);

		/* test for doneness */
		/* occasionally an iteration produces no improvement in chi-squ, but hasn't
			reached a true minimum, so make it keep going for a few interations of
			no improvement before actually stopping */
		if (( (fabs(fLastChi-fChi)/fChi < 0.00001 ) && (iNoImprov>1) )|| fChi <= fCaim || fabs(fLastChi-fChi) < 0.01) {
			/* no significant improvement... stop */
			break;
		}
		if (fabs(fLastChi-fChi)/fChi < 0.00001 || fabs(fLastChi-fChi) < 0.01){
			iNoImprov++;
		}
		else {
			iNoImprov = 0;
		}
		fLastChi = fChi;
	}

	sprintf(strMessage,"Finished after %d iterations. Last TEST value was: %g",iIteration,fTEST);
	TRACE(LOG_HIGH_PRI,strMessage);

	/* update final results */
	*pFinalChiSqu = fChi;
	*pFinalQ = fQ;
	*pFinalEntropy = fEnt;
	if (*pResultSource != NULL) {
		lv_free_image_struct(*pResultSource);
	}
	*pResultSource = lv_duplicate_image(g_pSrc);
	if (*pResultImage != NULL) {
		lv_free_image_struct(*pResultImage);
	}
	*pResultImage = lv_duplicate_image(pProjectedImage);

EXIT:
	/* free working space */
	for (i=0; i<NUM_SEARCH_DIR; i++) {
		if (pLinIndVector[i] != NULL) free(pLinIndVector[i]);
		if (pSearchVector[i] != NULL) free(pSearchVector[i]);
	}

	if (pGradGradC != NULL) free(pGradGradC);
	if (pGradS != NULL) free(pGradS);
	if (pGradC != NULL) free(pGradC);
	if (g_pTempImg != NULL) {
		lv_free_image_struct(g_pTempImg);
		g_pTempImg=NULL;
	}
	if (g_pSrc != NULL) {
		lv_free_image_struct(g_pSrc);
		g_pSrc=NULL;
	}
	if (g_pRevPSF != NULL) {
		lv_free_image_struct(g_pRevPSF);
		g_pRevPSF = NULL;
	}
	if (pPSF != NULL) lv_free_image_struct(pPSF);
    g_pVariance = NULL;

	TRACE_OUT;
	return iStatus;
}


/*************************
**************************/
static int	mem_CalcLinIndepModels(real_t *pSearchVectors[], real_t *pLinIndVectors[], real_t *pGradGradC, real_t pM_Eigenvales[], int iUseSearchVec[NUM_SEARCH_DIR]) {
	int	iResult=0,i;
	double	fMetric[NUM_SEARCH_DIR*NUM_SEARCH_DIR];
	double	fTransMetric[NUM_SEARCH_DIR*NUM_SEARCH_DIR];
	double	fCurv[NUM_SEARCH_DIR*NUM_SEARCH_DIR];
	double	fMaxEval=0, mag_e1=0,mag_e2=0,mag_e3=0;
	gsl_matrix_view m,trans_m;
	gsl_vector_view	col1,col2,col3;
	gsl_vector      *eval=NULL,*m_eval=NULL;
	gsl_matrix      *evec=NULL;
	gsl_matrix      *mTemp=NULL;	/* a temp copy of the raw metric which will
										be destroyed by gsl_eigen_symmv */
	gsl_eigen_symmv_workspace       *workspace=NULL;
	gsl_eigen_symm_workspace       *workspace_M=NULL;

	TRACE_IN(mem_CalcLinIndepModels);

	/* make working space */
	m = gsl_matrix_view_array(fMetric,NUM_SEARCH_DIR,NUM_SEARCH_DIR);
	trans_m = gsl_matrix_view_array(fTransMetric,NUM_SEARCH_DIR,NUM_SEARCH_DIR);
	eval = gsl_vector_alloc(NUM_SEARCH_DIR);
	evec = gsl_matrix_alloc(NUM_SEARCH_DIR,NUM_SEARCH_DIR);
	mTemp = gsl_matrix_alloc(NUM_SEARCH_DIR,NUM_SEARCH_DIR);
	workspace = gsl_eigen_symmv_alloc(NUM_SEARCH_DIR);
	workspace_M = gsl_eigen_symm_alloc(NUM_SEARCH_DIR);

	/* init */
	for (i=0; i< NUM_SEARCH_DIR*NUM_SEARCH_DIR; i++) fCurv[i]=0;
	for (i=0; i< NUM_SEARCH_DIR; i++) {
		lv_ZeroRealArray(pLinIndVectors[i],g_iSizeSrc);
		pM_Eigenvales[i] = 0;
	}

    /* normalise the vectors */
    mag_e1 = sqrt(img_VecScalarProd(pSearchVectors[0],pSearchVectors[0],g_iSizeSrc));
    mag_e2 = sqrt(img_VecScalarProd(pSearchVectors[1],pSearchVectors[1],g_iSizeSrc));
    mag_e3 = sqrt(img_VecScalarProd(pSearchVectors[2],pSearchVectors[2],g_iSizeSrc));

    /* protect against unlikely div by zero */
    if (mag_e1 == 0.0) mag_e1 = 1.0;
    if (mag_e2 == 0.0) mag_e2 = 1.0;
    if (mag_e3 == 0.0) mag_e3 = 1.0;
    img_VecMultByScalar(pSearchVectors[0],1./mag_e1,g_iSizeSrc);
    img_VecMultByScalar(pSearchVectors[1],1./mag_e2,g_iSizeSrc);
    img_VecMultByScalar(pSearchVectors[2],1./mag_e3,g_iSizeSrc);

	/* calculate g, the metric */
	/* e1 dot e1 */
	fMetric[0]                  = img_VecScalarProd(pSearchVectors[0],pSearchVectors[0],g_iSizeSrc);
	/* e2 dot e2 */
	fMetric[NUM_SEARCH_DIR+1]   = img_VecScalarProd(pSearchVectors[1],pSearchVectors[1],g_iSizeSrc);
	/* e3 dot e3 */
	fMetric[NUM_SEARCH_DIR*2+2] = img_VecScalarProd(pSearchVectors[2],pSearchVectors[2],g_iSizeSrc);

	/* e1 dot e2 */
	fMetric[1] = fMetric[NUM_SEARCH_DIR+0] =
      img_VecScalarProd(pSearchVectors[0],pSearchVectors[1],g_iSizeSrc);
	/* e1 dot e3 */
	fMetric[2] = fMetric[NUM_SEARCH_DIR*2+0] =
      img_VecScalarProd(pSearchVectors[0],pSearchVectors[2],g_iSizeSrc);
	/* e2 dot e3 */
	fMetric[NUM_SEARCH_DIR*2+1] = fMetric[NUM_SEARCH_DIR+2] =
      img_VecScalarProd(pSearchVectors[1],pSearchVectors[2],g_iSizeSrc);

	sprintf(strMessage,"metric row 1: %g,\t%g,\t%g",fMetric[0],fMetric[1],fMetric[2]);
	TRACE(LOG_MED_PRI,strMessage);
	sprintf(strMessage,"metric row 2: %g,\t%g,\t%g",fMetric[3],fMetric[4],fMetric[5]);
	TRACE(LOG_MED_PRI,strMessage);
	sprintf(strMessage,"metric row 3: %g,\t%g,\t%g",fMetric[6],fMetric[7],fMetric[8]);
	TRACE(LOG_MED_PRI,strMessage);

	/* make a copy of the metric, because it will be desroyed in the next step */
	gsl_matrix_memcpy(mTemp,&m.matrix);

	/* find e-values and e-vectors of g (destroys values input matrix)*/
	/* eigenvectors are in columns of evec */
	gsl_eigen_symmv(mTemp, eval, evec, workspace);

	/* don't sort by eigenvalue- need to retain index ordering for
       calculations outside this function. . */

	sprintf(strMessage,"Eigenvalues: %g,\t%g,\t%g",gsl_vector_get(eval,0),gsl_vector_get(eval,1),gsl_vector_get(eval,2));
	TRACE(LOG_MED_PRI,strMessage);

	col1 = gsl_matrix_column(evec,0);
	col2 = gsl_matrix_column(evec,1);
	col3 = gsl_matrix_column(evec,2);

	sprintf(strMessage,"Raw eigenvectors: [%g,%g,%g],\t[%g,%g,%g],\t[%g,%g,%g]",
				gsl_vector_get(&col1.vector,0),gsl_vector_get(&col1.vector,1),gsl_vector_get(&col1.vector,2),
				gsl_vector_get(&col2.vector,0),gsl_vector_get(&col2.vector,1),gsl_vector_get(&col2.vector,2),
				gsl_vector_get(&col3.vector,0),gsl_vector_get(&col3.vector,1),gsl_vector_get(&col3.vector,2));
	TRACE(LOG_MED_PRI,strMessage);

	/* need to protect against linearly dependent eigenvectors here.
		This shows up as relatively small eigenvalues. small is < 1e-5 smaller
		than biggest eval roughly */
	for(i=0; i< NUM_SEARCH_DIR; i++) {
		double ftemp=0;
		ftemp = gsl_vector_get(eval,i);
		if (ftemp > fMaxEval) {
			fMaxEval = ftemp;
		}
	}

	for(i=0; i<NUM_SEARCH_DIR; i++) {
        iUseSearchVec[i] = 1;
		if (fabs(gsl_vector_get(eval,i)/fMaxEval) < SMALL_EVAL) {
			sprintf(strMessage,"Discarding eigenval %g and corresponding eigenvector",gsl_vector_get(eval,i));
			TRACE(LOG_MED_PRI,strMessage);
			iUseSearchVec[i] = 0;
		}
	}

	/* scale columns of evector matrix (the evectors themselves) by sqrt(e-value) of transformed g 
		so that evec # g # transpose(evec) = I */
	/* scale eigenvectors so that the metric, when transformed below, will be cartesian */
    /*
	gsl_vector_scale(&col1.vector,1./sqrt(fabs(gsl_vector_get(eval,0))));
	gsl_vector_scale(&col2.vector,1./sqrt(fabs(gsl_vector_get(eval,1))));
	gsl_vector_scale(&col3.vector,1./sqrt(fabs(gsl_vector_get(eval,2))));
    */

	/* transform metric. This is transpose(E)#g#E where E is matrix of eigen vectors
		and g is metric and # is matrix multiply */
	/* this is actually redundant, but useful for checking that everything is working OK */
	cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,NUM_SEARCH_DIR,NUM_SEARCH_DIR,NUM_SEARCH_DIR,1.0,fMetric,NUM_SEARCH_DIR,evec->data,NUM_SEARCH_DIR,0.0,mTemp->data,NUM_SEARCH_DIR);

	/* mult intermediate by transposed matrix of eigenvectors */
	cblas_dgemm(CblasRowMajor,CblasTrans,CblasNoTrans,NUM_SEARCH_DIR,NUM_SEARCH_DIR,NUM_SEARCH_DIR,1.0,evec->data,NUM_SEARCH_DIR,mTemp->data,NUM_SEARCH_DIR,0.0,fTransMetric,NUM_SEARCH_DIR);

	sprintf(strMessage,"transformed metric row 1: %g,\t%g,\t%g",fTransMetric[0],fTransMetric[1],fTransMetric[2]);
	TRACE(LOG_MED_PRI,strMessage);
	sprintf(strMessage,"transformed metric row 2: %g,\t%g,\t%g",fTransMetric[3],fTransMetric[4],fTransMetric[5]);
	TRACE(LOG_MED_PRI,strMessage);
	sprintf(strMessage,"transformed metric row 3: %g,\t%g,\t%g",fTransMetric[6],fTransMetric[7],fTransMetric[8]);
	TRACE(LOG_MED_PRI,strMessage);

	/* calculate new, orthonormal search directions using scaled eigenvec matrix */
    for (i=0; i< g_iSizeSrc; i++) {
        if(iUseSearchVec[0]) {
            pLinIndVectors[0][i] = pSearchVectors[0][i]*gsl_vector_get(&col1.vector,0)
                                  +pSearchVectors[1][i]*gsl_vector_get(&col1.vector,1)
                                  +pSearchVectors[2][i]*gsl_vector_get(&col1.vector,2);
        }
        if(iUseSearchVec[1]) {
            pLinIndVectors[1][i] = pSearchVectors[0][i]*gsl_vector_get(&col2.vector,0)
                                  +pSearchVectors[1][i]*gsl_vector_get(&col2.vector,1)
                                  +pSearchVectors[2][i]*gsl_vector_get(&col2.vector,2);
        }
        if(iUseSearchVec[2]) {
            pLinIndVectors[2][i] = pSearchVectors[0][i]*gsl_vector_get(&col3.vector,0)
                                  +pSearchVectors[1][i]*gsl_vector_get(&col3.vector,1)
                                  +pSearchVectors[2][i]*gsl_vector_get(&col3.vector,2);
        }
    }


	/* calculate the "M" metric which represents the local curvature of the
		chi-squ surface. This is e_mu dot gradgradC dot e_nu*/
	for (i=0; i< g_iSizeSrc; i++) {
		fCurv[0]                            += pLinIndVectors[0][i]*pLinIndVectors[0][i]*pGradGradC[i];
		fCurv[1] = (fCurv[NUM_SEARCH_DIR]   += pLinIndVectors[0][i]*pLinIndVectors[1][i]*pGradGradC[i]);
		fCurv[2] = (fCurv[NUM_SEARCH_DIR*2] += pLinIndVectors[0][i]*pLinIndVectors[2][i]*pGradGradC[i]);
		fCurv[NUM_SEARCH_DIR+1]             += pLinIndVectors[1][i]*pLinIndVectors[1][i]*pGradGradC[i];
		fCurv[NUM_SEARCH_DIR+2] = (fCurv[NUM_SEARCH_DIR*2+1]
                                            += pLinIndVectors[1][i]*pLinIndVectors[2][i]*pGradGradC[i]);
		fCurv[NUM_SEARCH_DIR*2 + 2]         += pLinIndVectors[2][i]*pLinIndVectors[2][i]*pGradGradC[i];
	}

	sprintf(strMessage,"M row 1: %g, %g, %g",fCurv[0],fCurv[1],fCurv[2]);
	TRACE(LOG_MED_PRI,strMessage);
	sprintf(strMessage,"M row 2: %g, %g, %g",fCurv[3],fCurv[4],fCurv[5]);
	TRACE(LOG_MED_PRI,strMessage);
	sprintf(strMessage,"M row 3: %g, %g, %g",fCurv[6],fCurv[7],fCurv[8]);
	TRACE(LOG_MED_PRI,strMessage);

	/* now calculate the eigenvalues of the M metric. Destroys the
       matrix.*/
    /* this is redundant. M is already diagonal. Just pull the
       eigenvalues from the diagonal elements. Also, ordering of
       eigenvals/vecs does not correspond to anything in the original array
       m_eval = gsl_vector_alloc(NUM_SEARCH_DIR);
       m = gsl_matrix_view_array_with_tda(fCurv,NUM_SEARCH_DIR,NUM_SEARCH_DIR,NUM_SEARCH_DIR);
       gsl_eigen_symm(&m.matrix, m_eval, workspace_M);
    */
	for (i=0; i<NUM_SEARCH_DIR; i++) {
		pM_Eigenvales[i] = fCurv[NUM_SEARCH_DIR*i + i];
		sprintf(strMessage,"M eigenval %d: %g",i+1,pM_Eigenvales[i]);
		TRACE(LOG_MED_PRI,strMessage);
	}

/*EXIT:*/
	/* free vectors etc */
	if (mTemp != NULL) gsl_matrix_free(mTemp);
	if (eval != NULL) gsl_vector_free(eval);
	if (m_eval !=NULL) gsl_vector_free(m_eval);
	if (evec != NULL) gsl_matrix_free(evec);
	if (workspace!= NULL) gsl_eigen_symmv_free(workspace);
	if (workspace_M!= NULL) gsl_eigen_symm_free(workspace_M);

	TRACE_OUT;
	return iResult;
}


/*************************
**************************/
static	int	mem_CalcSearchDirs(lv_mapmatrix_t *pMap, lv_image_t *pData, lv_image_t *pImg, lv_image_t *pSrc, lv_image_t *pNoise, real_t *pSearchVector[NUM_SEARCH_DIR], real_t *pGradGradC, lv_image_t *pPSF,real_t *pGradS, real_t *pGradC, real_t *pTEST) {
	int	iResult=0,i;
	real_t	fMagS=0, fMagC=0;
	lv_image_t	*pTempGradC=NULL;

	TRACE_IN(mem_CalcSearchDirs);

	/* make space for gradC in image plane, which we then reverse project */
	pTempGradC = lv_duplicate_image(pData);
	if (pTempGradC==NULL) {
		iResult=-1;
		LOG_ERR("No malloc");
		goto EXIT;
	}

	lv_ZeroRealArray(pGradS,g_iSizeSrc);
	lv_ZeroRealArray(pGradC,g_iSizeSrc);

	/* calculate gradS */
	iResult = img_CalcEntropyDeriv(pSrc, pGradS);
	if (iResult != 0) {
		sprintf(strMessage,"img_CalcEntropyDeriv returned %d. Exiting...",iResult);
		goto EXIT;
	}

	/* calculate gradC */
	iResult = img_CalcChiSqDeriv(pData, pImg, pTempGradC->pImage, (g_imgCalcVarianceFlag == TRUE ? 0.0 : g_imgImgVariance), pNoise, pPSF);
	if (iResult != 0) {
		sprintf(strMessage,"img_CalcChiSqDeriv returned %d. Exiting...",iResult);
		goto EXIT;
	}

	if (g_iDebugImgs) img_DumpImage(pTempGradC->pImage,pTempGradC->pAxisSize,"gradCImg");

	/* reverse-project the gradient */
	iResult = img_ReverseProject(pMap, pTempGradC,pGradC);
	if (iResult !=0) {
		goto EXIT;
	}

	/* mult result by -1.0 because of difference in definition of
		C. Have been working with C = sum(D-F)^2/sig, but S&B
		definition is C=sum(F-D)^2/sig
	*/
	img_VecMultByScalar(pGradC,-1.0,g_iSizeSrc);
	img_VecMultByScalar(pGradS,-1.0,g_iSizeSrc);

	/* calculate the "magnitudes" of gradS and gradC */
	for (i=0; i< g_iSizeSrc; i++) {
		fMagS += pSrc->pImage[i]*pGradS[i]*pGradS[i];
		fMagC += pSrc->pImage[i]*pGradC[i]*pGradC[i];
	}
	fMagS = sqrt(fMagS);
	fMagC = sqrt(fMagC);

	if (fMagS == 0 || fMagC == 0) {
		sprintf(strMessage,"ERROR: fMagS (%g) or fMagC (%g) is zero. Can't calculate search directions",fMagS,fMagC);
		LOG_ERR(strMessage);
		iResult =1;
		goto EXIT;
	}

	/* calculate the search vectors */
	/* calculate (also) and print the "TEST" parameter which measures the 
		degree of parallelism between search directions */
	*pTEST=0;
	for (i=0; i< g_iSizeSrc; i++) {
		real_t	fTemp =0;
		pSearchVector[0][i] = pSrc->pImage[i]*pGradS[i];
		pSearchVector[1][i] = pSrc->pImage[i]*pGradC[i];
		pSearchVector[2][i] = pSrc->pImage[i]*pGradGradC[i]*(pSearchVector[0][i]/fMagS - pSearchVector[1][i]/fMagC);
		fTemp = (pGradS[i]/(fMagS) - pGradC[i]/fMagC);
		*pTEST += fTemp*fTemp;
	}

	sprintf(strMessage,"Mag S: %g. Mag C: %g.",fMagS,fMagC);
	TRACE(LOG_MED_PRI,strMessage);

EXIT:
	if(pTempGradC != NULL) lv_free_image_struct(pTempGradC);
	TRACE_OUT;
	return iResult;
}


/*************************
**************************/
static int mem_CalcGradGradC(real_t	*pGradGradC, lv_image_t *pNoise, lv_mapmatrix_t *pMap, lv_image_t *pPSF, lv_image_t *pRevPSF) {
	int	iResult=0, i;
	lv_image_t	*pTemp=NULL,*pTemp2=NULL;

	TRACE_IN(mem_CalcGradGradC);

	if (pNoise==NULL) {
		/* no noise file supplied. Just use the global varaince estimate */
		pTemp = lv_create_image_struct(pMap->dimension[0], pMap->dimension[1], IMGTYPE, g_PixelResn);
		if (pTemp==NULL) {
			LOG_ERR("No malloc");
			iResult =-1;
			goto EXIT;
		}
		for (i=0; i< g_iSizeImg; i++) {
			pTemp->pImage[i] = 2.0/g_imgImgVariance;
		}

	}
	else {
		real_t	fMean=0;
		int		iCount=0;

		pTemp = lv_duplicate_image(pNoise);
		pTemp2 = lv_duplicate_image(pNoise);
		if (pTemp == NULL || pTemp2==NULL) {
			LOG_ERR("No malloc")
			iResult =-1;
			goto EXIT;
		}

		/* fix variance value which indicate bad pixels */
		for (i=0; i< g_iSizeImg; i++) {
			if (pTemp->pImage[i] < 998.) {
				fMean+= pTemp->pImage[i];
				iCount++;
			}
		}
		fMean /=iCount;
		for (i=0; i< g_iSizeImg; i++) {
			if (pTemp->pImage[i]>998.) pTemp->pImage[i]=fMean;
		}

	/* experiment... 
		for (i=0; i< pPSF->pAxisSize[0]*pPSF->pAxisSize[1]; i++) {
			if (pPSF->pImage[i] > fMax) fMax = pPSF->pImage[i];
		}
		fMax = 1./fMax;
	*/

		/* convolve with PSF */
		if (g_iDebugImgs) img_DumpImage(pTemp->pImage,pTemp->pAxisSize,"gradGradCImg_pre");
		img_ConvolveImgsWithFFT(pTemp, pPSF, pTemp2);
		if (g_iDebugImgs) img_DumpImage(pTemp2->pImage,pTemp2->pAxisSize,"gradGradCImg_mid");
		img_ConvolveImgsWithFFT(pTemp2, pRevPSF, pTemp);
		if (g_iDebugImgs) img_DumpImage(pTemp->pImage,pTemp->pAxisSize,"gradGradCImg");

		/* now transform the PSF-convolved image into the 2.0/variance 
			 which is required for grad(gradC) */
		for (i=0; i< g_iSizeImg; i++) {
			if (pTemp->pImage[i] <= 0) {
				LOG_ERR("ERROR: negative/zero values in variance file");
				iResult = 1;
				goto EXIT;
			}
			/*
			pTemp->pImage[i] = 2.0/(pTemp->pImage[i]*fMax);
			*/
			pTemp->pImage[i] = 2.0/(pTemp->pImage[i]);
		}
	}

	/* reverse project the 2/noise pixel values */
	img_ReverseProject(pMap,pTemp,pGradGradC);

EXIT:
	if (pTemp!=NULL) lv_free_image_struct(pTemp);
	if (pTemp2!=NULL) lv_free_image_struct(pTemp2);
	TRACE_OUT;
	return iResult;
}

/*************************
**************************/
static int	mem_CalcQ(lv_mapmatrix_t *pMap,lv_image_t *pSrc, lv_image_t *pImg, real_t *pResultChiSqu, real_t *pResEntropy, lv_image_t *pPSF){
	int	iStatus=0;

	TRACE_IN(mem_CalcQ);

	/* calc the entropy */
	img_CalcImgEntropy(pSrc,pResEntropy,pSrc->fDefaultVal);

	/* forward-project the source then calc the chi-squ */
	lv_ZeroRealArray(g_pTempImg->pImage,img_CalcImgSize(g_pTempImg));
	lv_projectSourceThruMapMatrix(pSrc, g_pTempImg, pMap);
	if (pPSF != NULL) {
		sprintf(strMessage,"Convolving mapped image with PSF.");
		TRACE(LOG_LOW_PRI,strMessage);
		if (g_iDebugImgs) img_DumpImage(g_pTempImg->pImage,g_pTempImg->pAxisSize,"imgIter_pre");
        img_ConvolveImgsWithFFT(g_pTempImg,pPSF,pImg);
		if (g_iDebugImgs) img_DumpImage(pImg->pImage,pImg->pAxisSize,"imgIter_post");
    }

	if (g_bUseMultImgPixOnly) {
		*pResultChiSqu = img_CalcChiSquared(g_pData,pImg, (g_imgCalcVarianceFlag ? 0.0 :g_imgImgVariance), g_pVariance,pMap->pMultImgPix);
	}
	else {
		*pResultChiSqu = img_CalcChiSquared(g_pData,pImg, (g_imgCalcVarianceFlag ? 0.0 :g_imgImgVariance), g_pVariance,pImg->pMask);
	}
	
/*EXIT:*/
	TRACE_OUT;
	return iStatus;
}

