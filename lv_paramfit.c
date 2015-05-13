/*************************
Code to implement the downhill simplex method of parameter
optimisation as part of the Lensview grav lensing software.
Copyright (C) 2006. Randall Wayth.
$Id: lv_paramfit.c,v 1.5 2008/10/31 20:43:10 rwayth Exp rwayth $
**************************/
#include	<unistd.h>
#include	<stdio.h>
#include	<stdlib.h>
#include	<string.h>
#include	<math.h>
#include	"gsl/gsl_multimin.h"
#include	"common.h"
#include	"log.h"
#include	"lv_common.h"
#include	"lv_mem.h"
#include	"lv_paramfit.h"
#include	"lv_lens.h"
#include    "lv_image.h"

#define		SMALL_NUM	1.e-10
#define		TOLERANCE	1e-3
#define		MAX_ITERATIONS	200

/* structure to hold the params required for minimisation */
typedef struct _minimiser_params {
    lv_image_t *g_pPSF, *g_RevPSF, *g_Variance;
    lv_lensmodel_t *pLens;
    lv_mapmatrix_t  *g_pWeightMatrix;
    lv_image_t	*g_pSourceOriginal, *g_pMappedImage, *g_pData;
    int		g_iMaxIterations,iNParams;
    real_t *pResChiSqu, *pResEntropy, g_fResQ;
    lv_image_t **pBestSrc, **pBestImg, *g_pNoiseImg;
    float g_fSrcDefaultVal,fCentreSearchRange;
    bool g_SearchCentre;
} minimiser_params;


static	double lp_doMEM_wrap(const gsl_vector *v, void *voidparams);
static	int lp_gsl_min_wrap(minimiser_params *params);
static  int countIterations(lv_lensmodel_t *pLens);
static  int updatePararms(lv_lensmodel_t *pLens);


/***************************
Function:   lp_MinFinder
Description:    Entry point for amoeba style minimiser.
Arguments:
Returns:    0: success.
****************************/
int	lp_MinFinder(lv_lensmodel_t *pLens, lv_image_t *pSource, lv_image_t *pMappedImage, lv_image_t *pData,
                int iMaxIterations,real_t *pResChiSqu, real_t *pResEntropy, lv_image_t **pBestSrc, lv_image_t **pBestImg,
                lv_image_t *pPsfImg, lv_image_t *pNoiseImg, int bSearchCentre, float fSrcDefaultVal) {
                
    minimiser_params params;
    /* initialise the minimiser parameters */
    
	params.g_pPSF = pPsfImg;
	params.g_RevPSF=NULL;
	params.g_pWeightMatrix = NULL;
	params.g_Variance =  pNoiseImg;
	params.pLens = pLens;
	params.g_pSourceOriginal = pSource;
	params.g_pMappedImage = pMappedImage;
	params.g_pData = pData;
	params.g_iMaxIterations = iMaxIterations;
	params.pResChiSqu = pResChiSqu;
	params.pResEntropy = pResEntropy;
	params.pBestSrc = pBestSrc;
	params.pBestImg = pBestImg;
	params.g_pNoiseImg = pNoiseImg;
	params.g_SearchCentre = bSearchCentre;
	params.g_fSrcDefaultVal = fSrcDefaultVal;
	params.fCentreSearchRange = 0.25;
	params.iNParams = 0;

	return lp_gsl_min_wrap(&params);
}


/***************************
Function:   lp_gsl_min_wrap
Description:    Wrapper function that sets up the GSL amoeba style minimiser then calls the minimiser.
Arguments:
Returns:    0: success.
****************************/
static int	lp_gsl_min_wrap(minimiser_params *params) {
	size_t iter = 0;
	int status=0,i,j,iParam=0,iCount=0,iOffset=0,k=0;
	const gsl_multimin_fminimizer_type *T;
	gsl_multimin_fminimizer *s=NULL;
	gsl_vector *x=NULL,*step=NULL;
	gsl_multimin_function my_func;
	double	fSize=0,fTemp=0;

	TRACE_IN(lp_gsl_min_wrap);

	/* how many params are there? */
	for (i=0; i< params->pLens->iNumComponents; i++) {
		params->iNParams += params->pLens->pComponent[i].iNumParams;
	}
	if (params->g_SearchCentre) params->iNParams += 2;

	sprintf(strMessage,"There are %d parameters. Searching centre: %d",params->iNParams,(int)params->g_SearchCentre);
	TRACE(LOG_HIGH_PRI,strMessage);

	my_func.f = &lp_doMEM_wrap;
	my_func.n = params->iNParams;
	my_func.params = params;

	/* Set up initial guess. Just take the middle of the param range */
	x = gsl_vector_alloc(params->iNParams);
	iParam=0;
	if (params->g_SearchCentre) {
		iParam=2;
		gsl_vector_set(x, 0,params->pLens->pComponent[0].fXoffset);
		gsl_vector_set(x, 1,params->pLens->pComponent[0].fYoffset);
	}
	for (j=0; j< params->pLens->iNumComponents; j++) {
		for (i=0; i< params->pLens->pComponent[j].iNumParams ; i++) {
			fTemp = (params->pLens->pComponent[j].fParamto[i] + params->pLens->pComponent[j].fParamfrom[i])*0.5;
			gsl_vector_set (x, iParam,fTemp);
			sprintf(strMessage,"Param %d. Staring at: %g",iParam,fTemp);
			TRACE(LOG_HIGH_PRI,strMessage);
			iParam++;
		}
	}

	/* set step sizes */
	step = gsl_vector_alloc(params->iNParams);
	iParam=0;
	if (params->g_SearchCentre) {
		iParam=2;
		gsl_vector_set(step, 0,params->fCentreSearchRange);
		gsl_vector_set(step, 1,params->fCentreSearchRange);
	}
	for (j=0; j< params->pLens->iNumComponents; j++) {
		for (i=0; i< params->pLens->pComponent[j].iNumParams ; i++) {
			fTemp = (params->pLens->pComponent[j].fParamto[i] - params->pLens->pComponent[j].fParamfrom[i])*0.5;
			gsl_vector_set (step, iParam,fTemp);
			sprintf(strMessage,"Param %d. inital step: %g",iParam,fTemp);
			TRACE(LOG_HIGH_PRI,strMessage);
			iParam++;
		}
	}

	T = gsl_multimin_fminimizer_nmsimplex;
	s = gsl_multimin_fminimizer_alloc(T, params->iNParams);

	gsl_multimin_fminimizer_set(s, &my_func, x, step);

	do {
		iter++;
		status = gsl_multimin_fminimizer_iterate(s);

		if (status) break;

		fSize = gsl_multimin_fminimizer_size (s);
		status = gsl_multimin_test_size (fSize, 1e-2);
		sprintf(strMessage,"Finished iteration %d. chisqu: %g. Simplex size: %g",(int)iter,s->fval,fSize);
		TRACE(LOG_HIGH_PRI,strMessage);

	} while (status == GSL_CONTINUE && iter < 1000);

	if (status == GSL_SUCCESS) {
		iCount = sprintf(strMessage,"Finished minimisation. chisqu: %g entropy: %g Best fit: lensoffset x: %g offsety: %g",
		                s->fval,(double)*(params->pResEntropy),params->pLens->pComponent[0].fXoffset,params->pLens->pComponent[0].fYoffset);

		iOffset += iCount;
		if (params->g_SearchCentre) k=2;
		for (j=0; j< params->pLens->iNumComponents; j++) {
			for (i=0; i< params->pLens->pComponent[j].iNumParams; i++) {
				iCount = sprintf(strMessage+iOffset," %s: %g",params->pLens->pComponent[j].strParamName[i],gsl_vector_get(s->x,k++));
				iOffset += iCount;
			}
		}
		TRACE(LOG_HIGH_PRI,strMessage);
	}

//EXIT:
	if (s!=NULL) gsl_multimin_fminimizer_free(s);
	if (x!=NULL) gsl_vector_free(x);
	if (step!=NULL) gsl_vector_free(step);
	TRACE_OUT;

	return 0;
}


/***************************
Function:   lp_doMEM_wrap
Description: This function is called by GSL code that runs the minimiser. It unpacks lens params, creates mapping matrix
             and calls MEM function for the lens model. The resulting chisqu is returned to the minimiser.
Arguments:
Returns:    0: success.
****************************/
static double lp_doMEM_wrap(const gsl_vector *v, void *voidparams) {
	int iStatus=0, iParam=0,iNumImgPix=0, iNumSrcPix=0,i,j;
	int	iOffset=0, iCount=0;
	real_t fReducedChiSqu=0;
	lv_image_t *pSource=NULL;
    minimiser_params *params;
    
	TRACE_IN(lp_doMEM_wrap);
	
	params = voidparams;
	/* extract the lens model params from the vector */
	iParam=0;
	if (params->g_SearchCentre) {
		iParam=2;
		params->pLens->pComponent[0].fXoffset = gsl_vector_get(v,0);
		params->pLens->pComponent[0].fYoffset = gsl_vector_get(v,1);
	}
	for (j=0; j< params->pLens->iNumComponents; j++) {
		for (i=0; i< params->pLens->pComponent[j].iNumParams ; i++) {
			params->pLens->pComponent[j].fParameter[i] = gsl_vector_get(v,iParam);
			iParam++;
		}
		/* by default the external shear is at the grid centre.
		 * need to make it for the lens centre */
		if (params->pLens->pComponent[j].iType == LM_EXTSHEAR) {
			params->pLens->pComponent[j].fXoffset = params->pLens->pComponent[0].fXoffset;
			params->pLens->pComponent[j].fYoffset = params->pLens->pComponent[0].fYoffset;
		}
	}

	iCount = sprintf(strMessage,"Starting minimisation with lensoffset x: %g, offsety: %g",
	                params->pLens->pComponent[0].fXoffset,params->pLens->pComponent[0].fYoffset);
	iOffset += iCount;

	for (j=0; j< params->pLens->iNumComponents; j++) {
		for (i=0; i< params->pLens->pComponent[j].iNumParams; i++) {
			iCount = sprintf(strMessage+iOffset," %s: %g",params->pLens->pComponent[j].strParamName[i],params->pLens->pComponent[j].fParameter[i]);
			iOffset += iCount;
		}
	}
	TRACE(LOG_HIGH_PRI,strMessage);

	
	/* free previous mapping matrix ? */
	if(params->g_pWeightMatrix != NULL) lv_freeMappingMatrix(params->g_pWeightMatrix);

	/* copy the original source */
	pSource = lv_duplicate_image(params->g_pSourceOriginal);

	/* set up lens mapping matrix */
	params->g_pWeightMatrix = lv_allocMappingMatrix(params->g_pMappedImage->pAxisSize[0],params->g_pMappedImage->pAxisSize[1],
	                                    pSource->pAxisSize[0],pSource->pAxisSize[1]);
	if (params->g_pWeightMatrix == NULL) {
		fprintf(stderr,"lv_allocMappingMatrix failed.\n");
		LOG_ERR("lv_allocMappingMatrix failed. Exiting program...");
		exit(iStatus);
	}
	
	iStatus = lv_createMappingMatrix(params->pLens, params->g_pWeightMatrix, pSource, params->g_pMappedImage);
	if (iStatus != 0) {
		LOG_ERR("lv_createMappingMatrix failed. Failing for this param combination.");
	}

	/* do the fit */
    if (iStatus ==0) {
	    iStatus = mem_DoMEM(params->pLens, params->g_pWeightMatrix, pSource, params->g_pMappedImage, params->g_pData,
	                    params->g_iMaxIterations, params->pResChiSqu, params->pResEntropy, params->pBestSrc,
	                    params->pBestImg, params->g_pPSF, params->g_pNoiseImg,params->g_fSrcDefaultVal,&(params->g_fResQ));
    }

	if (iStatus == 0) {
		sprintf(strMessage,"MaxEn returned %d Chi Squ: %g Entropy: %g Reduced chisqu: %g Num img pix: %d Num src pix: %d",
		                iStatus,(double) *(params->pResChiSqu),(double) *(params->pResEntropy),(double)fReducedChiSqu,
		                iNumImgPix,iNumSrcPix);
		TRACE(LOG_HIGH_PRI,strMessage);
	}
	else {
		sprintf(strMessage,"WARNING. MaxEn returned %d. setting result to 1e9 and continuing.",iStatus);
		*(params->pResChiSqu) = 1e9;
		LOG_ERR(strMessage);
	}
//EXIT:
	TRACE_OUT;
	if (pSource != NULL) lv_free_image_struct(pSource);
	return *(params->pResChiSqu);
}


/***************************
Function:   lp_ParamSweep
Description: Wrapper function to perform a parameter sweep and calc the chisqu for each set of lens model params.
Arguments:
Returns:    0: success.
****************************/
int lp_ParamSweep(lv_lensmodel_t *pLens, lv_image_t *pSourceOriginal, lv_image_t *pMappedImage, lv_image_t *pData,
                int iMaxIterations,real_t *pBestChiSqu, real_t *pBestEntropy, lv_image_t **pBestSrc,
                lv_image_t **pBestImg, lv_image_t *pPsfImg, lv_image_t *pNoiseImg, int bSearchCentre,
                float fSearchCentreRange, float fSearchCentreStep, float fSrcDefaultVal,
                int iMethod, int iDumpImgs) {
    int iStatus=0,iIter=0,i,j;
    float	xcent=0,ycent=0,x_cent_offset=0,y_cent_offset=0,temp;
	lv_mapmatrix_t	*pWeightMatrix = NULL;
	lv_image_t *pSource=NULL,*pSrcResult=NULL,*pImgResult=NULL;

    TRACE_IN(lp_doParamSweep);
    
	/* keep an original copy of the soure */
	pSource = lv_duplicate_image(pSourceOriginal);
    *pBestChiSqu = 1e6;
    temp = img_CalcVariance(pData->pImage,img_CalcImgSize(pData));
    g_DataMean = img_CalcArrayTotal(pData->pImage,img_CalcImgSize(pData))/img_CalcImgSize(pData);
    sprintf(strMessage,"Data variance is %g. Mean is %g",temp,g_DataMean);
    TRACE(LOG_HIGH_PRI,strMessage);
    if (g_imgImgVariance == 0 ) g_imgImgVariance = temp;
    if (g_imgImgVariance ==0 && !g_imgCalcVarianceFlag) {
        LOG_ERR("ERROR: Global variance value is zero. Cannot continue. Exiting");
        iStatus=1;
        goto EXIT;
    }

    xcent = pLens->pComponent[0].fXoffset;
    ycent = pLens->pComponent[0].fYoffset;

    sprintf(strMessage,"Total expected paramter combinations: %d",countIterations(pLens));
    TRACE(LOG_HIGH_PRI,strMessage);

    if (bSearchCentre != FALSE) {
        if (fSearchCentreRange==0 || fSearchCentreStep==0) {
            sprintf(strMessage,"WARNING: search centre flag on, but range(%g)/step(%g) is zero. not searching centre",
				                fSearchCentreRange,fSearchCentreStep);
            LOG_ERR(strMessage);
            fSearchCentreRange=0;
            fSearchCentreStep=1;
        }
        else {
            int nstep=0;
            nstep = 2.0*fSearchCentreRange/fSearchCentreStep;
            sprintf(strMessage,"Cent search range: %g, step size: %g. Performing %d centre search steps for each lens model",
				        fSearchCentreRange,fSearchCentreStep,nstep);
            TRACE(LOG_HIGH_PRI,strMessage);
        }
    }
    else fSearchCentreStep=1.0;

    for (x_cent_offset=-fSearchCentreRange; x_cent_offset <= fSearchCentreRange; x_cent_offset+=fSearchCentreStep) {
    for (y_cent_offset=-fSearchCentreRange; y_cent_offset <= fSearchCentreRange; y_cent_offset+=fSearchCentreStep) {

    do {
        int	iOffset=0, iCount=0, iDOF=0, iNumImgPix=0, iNumSrcPix=0, iSizeSrc=0;
        real_t	fReducedChiSqu=0,fResChiSqu,fResEntropy,fResQ;

        iIter++;

		pLens->pComponent[0].fXoffset = xcent+x_cent_offset;
		pLens->pComponent[0].fYoffset = ycent+y_cent_offset;
		iCount = sprintf(strMessage,"Starting minimisation with lensoffset x: %g offsety: %g",
		                pLens->pComponent[0].fXoffset,pLens->pComponent[0].fYoffset);
		iOffset += iCount;

		for (j=0; j< pLens->iNumComponents; j++) {
			for (i=0; i< pLens->pComponent[j].iNumParams; i++) {
				iCount = sprintf(strMessage+iOffset," %s: %g",pLens->pComponent[j].strParamName[i],
				                pLens->pComponent[j].fParameter[i]);
				iOffset += iCount;
			}
		}

		TRACE(LOG_HIGH_PRI,strMessage);

		/* reset the mapping matrix for the next run with different parameters */
		if (pWeightMatrix != NULL) lv_freeMappingMatrix(pWeightMatrix);
		pWeightMatrix = lv_allocMappingMatrix(pMappedImage->pAxisSize[0],pMappedImage->pAxisSize[1],pSource->pAxisSize[0],
		                                        pSource->pAxisSize[1]);
		if (pWeightMatrix == NULL) {
			fprintf(stderr,"lv_allocMappingMatrix failed.\n");
			goto EXIT;
		}

		iStatus = lv_createMappingMatrix(pLens, pWeightMatrix, pSource, pMappedImage);
		if (iStatus != 0) {
			sprintf(strMessage,"lv_createMappingMatrix failed. invalid parameter val? Setting this param combination to large chi-squ...\n");
		    TRACE(LOG_ERROR_PRI,strMessage);
            fResChiSqu = 1e10;
			//goto EXIT;
            continue;
		}

		fResChiSqu = fResEntropy = fResQ = 0.0;

		/* note: MUST call img_DoConjGradient/mem_DoMEM with pMappedImage otherwise img plane pixel mask isn't defined */
		if (iMethod != 0) {
			iStatus = img_DoConjGradient(pLens, pWeightMatrix, pSource, pMappedImage, pData, iMaxIterations, &fResChiSqu,
			                             &fResEntropy, &pSrcResult, &pImgResult, pPsfImg, pNoiseImg);
		}
		else {
			iStatus = mem_DoMEM(pLens, pWeightMatrix, pSource, pMappedImage, pData, iMaxIterations, &fResChiSqu,
			                    &fResEntropy, &pSrcResult, &pImgResult, pPsfImg, pNoiseImg,fSrcDefaultVal,&fResQ);
		}

		if (iStatus == 0) {
			iSizeSrc = img_CalcImgSize(pSrcResult);
			iNumImgPix = img_CountActiveImgPixels(pMappedImage);
			iNumSrcPix = img_CountNonZeroElements(pSrcResult->pImage,iSizeSrc,fSrcDefaultVal);

			iDOF = iNumImgPix-iNumSrcPix - pLens->iNumParameters;
			fReducedChiSqu = fResChiSqu/(iDOF==0 ? 1 : iDOF);

			sprintf(strMessage,"MaxEn returned %d Chi Squ: %g Entropy: %g Reduced chisqu: %g Num img pix: %d Num src pix: %d",
			                iStatus,fResChiSqu,fResEntropy,fReducedChiSqu,iNumImgPix,iNumSrcPix);
			TRACE(LOG_HIGH_PRI,strMessage);
			
			/* keep track of the best result */
            if (*pBestChiSqu > fResChiSqu) {
				*pBestChiSqu = fResChiSqu;
				*pBestEntropy = fResEntropy;
				sprintf(strMessage,"Best Chisqu of %g. ent: %g", fResChiSqu,fResEntropy);
				TRACE(LOG_HIGH_PRI,strMessage);

				if (*pBestSrc != NULL) {
				    lv_free_image_struct(*pBestSrc);
				}
				*pBestSrc = lv_duplicate_image(pSrcResult);
				if (*pBestImg != NULL) {
				    lv_free_image_struct(*pBestImg);
				}
				*pBestImg = lv_duplicate_image(pImgResult);
			}

			/* dump the source/image if it is within 10 % of the best chi */
			if (iDumpImgs && ((fResChiSqu - *pBestChiSqu)/(*pBestChiSqu) ) < 0.1 ) {
			    char strOutputImageFile[FILENAME_MAX];
				sprintf(strMessage,"Dumping images for iteration %d",iIter);
				TRACE(LOG_HIGH_PRI,strMessage);
				sprintf(strOutputImageFile,"sourcedump_main%d.fits",iIter);
				lv_write_image_to_file(pSrcResult,strOutputImageFile,TRUE);
				sprintf(strOutputImageFile,"imagedump_main%d.fits",iIter);
				lv_write_image_to_file(pImgResult,strOutputImageFile,TRUE);
			}
		}
		else {
			sprintf(strMessage,"fitting process returned %d.",iStatus);
			TRACE(LOG_HIGH_PRI,strMessage);
		}

		/* all done with this set of params. reset the source */
		lv_free_image_struct(pSource);
		pSource = lv_duplicate_image(pSourceOriginal);
		lv_free_image_struct(pSrcResult); pSrcResult=NULL;
		lv_free_image_struct(pImgResult); pImgResult=NULL;
	} while (updatePararms(pLens) == 0) ;

	}
	}
EXIT:
    if(pSource !=NULL) lv_free_image_struct(pSource);
	if (pWeightMatrix != NULL) lv_freeMappingMatrix(pWeightMatrix);
	TRACE_OUT;
    return iStatus;
}


/******************************
 *****************************/
int	countIterations(lv_lensmodel_t *pLens) {
	int	count=1,i,j;
	for (i=0; i< pLens->iNumComponents; i++) {
		for (j=0; j< pLens->pComponent[i].iNumParams; j++) {
			if (pLens->pComponent[i].fParaminc[j] != 0) {
				count *= (1+ ((pLens->pComponent[i].fParamto[j]-pLens->pComponent[i].fParamfrom[j])/pLens->pComponent[i].fParaminc[j]));
/*
*				printf("%g,%g,%g means %d iterations for %s\n",
*								pLens->pComponent[i].fParamfrom[j],
*								pLens->pComponent[i].fParamto[j],
*								pLens->pComponent[i].fParaminc[j],
*								iter,pLens->pComponent[i].strParamName[j]);
*/
			}

		}
	}
	return count;
}


/* update the lens model parameters. If a parameter has reached its upper limit, then reset it
 * and update the next one down... and so on */
int	updatePararms(lv_lensmodel_t *pLens) {
	int i, iCurrParam, iDone=1,iNextParam=1;

	for(i = pLens->iNumComponents-1; i>=0; i--) {
		iCurrParam = pLens->pComponent[i].iNumParams-1;
/*
*		printf("Updating current param: %d of %d\n",iCurrParam,pLens->pComponent[i].iNumParams-1);
*/
		iNextParam=1;
		while(iNextParam) {
			pLens->pComponent[i].fParameter[iCurrParam] += pLens->pComponent[i].fParaminc[iCurrParam];
			if (pLens->pComponent[i].fParameter[iCurrParam] > (pLens->pComponent[i].fParamto[iCurrParam])*(1.0+1e-5) ||
				pLens->pComponent[i].fParamfrom[iCurrParam]==pLens->pComponent[i].fParamto[iCurrParam] ) {

				pLens->pComponent[i].fParameter[iCurrParam] = pLens->pComponent[i].fParamfrom[iCurrParam];
/*
*				printf("Resetting param %d to start\n",iCurrParam);
*/
				iCurrParam--;
				if (iCurrParam < 0) {
					iNextParam =0;
				}
			}
			else {
/*
*				printf("Updated param %d.\n",iCurrParam);
*/
				iDone = 0;
				goto EXIT;
			}
		}
	}
EXIT:
	return iDone;
}


