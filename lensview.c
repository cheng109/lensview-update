/* lensview.c: entry point and main control program for the
	"lensview" grav lens modelling software
Copyright (C) 2006. Randall Wayth.
$Id: lensview.c,v 1.20 2008/10/31 21:05:49 rwayth Exp rwayth $
*/
#include	<stdio.h>
#include	<stdlib.h>
#include	<time.h>
#include	<string.h>
#include	<math.h>
#include	<unistd.h>
#include	<sys/utsname.h>
#include	<limits.h>
#include	"fitsio.h"
#include	"common.h"
#include	"log.h"
#include	"lv_common.h"
#include	"lv_lens.h"
#include	"parseopts.h"
#include	"lv_image.h"
#include	"lv_mem.h"
#include	"lv_critline.h"
#include	"lv_paramfit.h"

/****************************
Public Global Variables
*****************************/
lv_axissize_t	g_iSzImgx=0, g_iSzImgy=0;

void printUsage(int argc, char * const argv[]);

/**************************
Private Global Variables
**************************/
//extern char    strParamFile[PATH_MAX] = "comps.txt";


int	main(int argc, char * const argv[]) {

	char	strOutputMappedFileName[PATH_MAX] = "model_img.fits", strOutputImageFile[PATH_MAX]="test_image.fits";
	char	strOutputSourceFile[PATH_MAX]="model_src.fits", strOutSrcMagFilename[] = "mag_src.fits";
	char	*strOptionString,strTempLogFilePath[PATH_MAX]="./", strSourceFileName[PATH_MAX] = "";
	char	strPsfFilename[PATH_MAX] = "", strDataFileName[PATH_MAX] = "", strNoiseFileName[PATH_MAX]="";
	char	strOutImgMagFilename[] = "mag_img_inv.fits";
    char    strParamFile[PATH_MAX] = "comps.txt";
	char	strMaskName[PATH_MAX]="";
	int		iStatus=0,iRayTraceOnly=FALSE,iNicePri=0,iMinimise=FALSE, iDumpImgs = FALSE;
	int		iMethod=0;
	int		iMaxIterations = 100, iMakeSrcInvMagMap=0;
	int		iNormaliseMax=0, iSearchCentre=0, iUseMinFinder=FALSE;
	float	x_axisoffset =0,y_axisoffset=0;
	float	fPixelResnRatio =1.0,fSrcDefaultVal=0,fConstVar=0;
	real_t	fResEntropy = 0.0, fResChiSqu = 0.0;
	float	temp=0.0,fSearchCentreRange=0,fSearchCentreStep=0;
	lv_image_t	*pSourceOriginal = NULL, *pSource = NULL, *pProjImage=NULL, *pMappedImage=NULL, *pPsfImg=NULL;
	lv_image_t	*pReadImg=NULL, *pNoiseImg=NULL, *pTemp=NULL;
	lv_image_t	*pBestSrc = NULL, *pBestImg = NULL, *pInvMag=NULL, *pMask=NULL;
	lv_lensmodel_t	*pLens = NULL;
	lv_mapmatrix_t	*pWeightMatrix = NULL;
	struct	utsname	mach_name;
	cl_pix_list cl_list_head = {0,0,0,0,NULL}, caust_list_head = {0,0,0,0,NULL};
	int	i; 
    int numJobs = 1;
/*
	printf("Hello world...\n");
	printf("sizeof int: %d, sizeof long: %d, sizeof short: %d, sizeof *: %d, float: %d, double: %d\n",(int)sizeof(int), (int)sizeof(long),(int) sizeof(short),(int) sizeof(void *),(int) sizeof(float), (int)sizeof(double));
	exit(0);
*/

	g_iTracePriority = LOG_LOW_PRI;

	strOptionString = "%-imgfile%s%-logfilepath%s%-srcxoffset%f%-srcyoffset%f%-tracelevel%i%-raytraceonly%t%-sourcefile%s%-psffile%s%-datafile%s%-dofit%t%-nice%i%-pixelratio%f%-maxiter%i%-pixelres%f%-noisefile%s%-makemag%t%-dumpimg%t%-paramfile%s%-fixlambda%f%-srcdefaultval%f%-mask%s%-useconjgrad%t%-normalisepsfmax%t%-debugimgs%t%-targetchisqu%f%-usemultimgpix%t%-constvariance%f%-searchcentre%t%-searchcentrerange%f%-searchcentrestep%f%-useminfinder%t%-numJobs%d%t";

	iStatus = parse_options(argc, argv,strOptionString ,strOutputImageFile, strTempLogFilePath, &x_axisoffset,&y_axisoffset,&g_iTracePriority,&iRayTraceOnly, strSourceFileName, strPsfFilename,strDataFileName,&iMinimise,&iNicePri,&fPixelResnRatio,&iMaxIterations, &g_PixelResn,strNoiseFileName,&iMakeSrcInvMagMap,&iDumpImgs,strParamFile,&g_FixedLambda,&fSrcDefaultVal,strMaskName,&iMethod,&iNormaliseMax,&g_iDebugImgs,&g_fTargetChiSqu,&g_bUseMultImgPixOnly,&fConstVar,&iSearchCentre,&fSearchCentreRange,&fSearchCentreStep,&iUseMinFinder);
    
    
    
    
    
    // Jun Cheng's update begin:
    char dynamicOutputModelImage[PATH_MAX] = "model_img_";
    strcat(dynamicOutputModelImage, strParamFile);
    strcat(dynamicOutputModelImage, ".fits");
   
    char dynamicOutputMaskImage[PATH_MAX] = "imgmask_";
    strcat(dynamicOutputMaskImage, strParamFile);
    strcat(dynamicOutputMaskImage, ".fits");
    
    char dynamicOutputSourceImage[PATH_MAX] = "model_src_";
    strcat(dynamicOutputSourceImage, strParamFile);
    strcat(dynamicOutputSourceImage, ".fits");
    
    char dynamicOutputResidulImage[PATH_MAX] = "model_res_";
    strcat(dynamicOutputResidulImage, strParamFile);
    strcat(dynamicOutputResidulImage, ".fits");

    // Jun Cheng's update end:
	
    
    if (argc == 1 || iStatus != 0 ) {
	  printUsage(argc,argv);
	}

	if (strcmp(strTempLogFilePath,"") != 0) {
		g_strLogFilePath = strdup(strTempLogFilePath);
	}
    if (strcmp(strParamFile,"") != 0) {
        g_paramFileName = strdup(strParamFile);
    }

	TRACE_IN(main);

	uname(&mach_name);
	sprintf(strMessage,"%s starting on machine %s",argv[0],mach_name.nodename);
	TRACE(LOG_HIGH_PRI,strMessage);

	/* be nice to other users */
	if (iNicePri > 0) {
		sprintf(strMessage,"Setting nice pri to %d",iNicePri);
		TRACE(LOG_HIGH_PRI,strMessage);
		nice(iNicePri);
	}

	if (g_FixedLambda != -1) {
		sprintf(strMessage,"Using fixed value for lambda: %g",g_FixedLambda);
		TRACE(LOG_HIGH_PRI,strMessage);
	}
	if (g_fTargetChiSqu >=0) {
		sprintf(strMessage,"Using user-defined target chi-squ: %g",g_fTargetChiSqu);
		TRACE(LOG_HIGH_PRI,strMessage);
	}

	if (g_PixelResn < 0.001) {
		LOG_ERR("WARNING. pixel resolution seems very small...");
	}

    /* read the data image */
	if (strDataFileName[0] != '\0') {
		sprintf(strMessage,"Reading data file <%s>",strDataFileName);
		TRACE(LOG_HIGH_PRI,strMessage);
		pReadImg = lv_read_img_file(strDataFileName,g_PixelResn);
		if (pReadImg == NULL) {
			sprintf(strMessage,"Failed to read FITS file %s.",strDataFileName);
			LOG_ERR(strMessage);
			fprintf(stderr,"%s\n",strMessage);
			exit(1);
		}
		g_iSzImgx = pReadImg->pAxisSize[0];
		g_iSzImgy = pReadImg->pAxisSize[1];

		/* read any external image plane mask */
		if(strMaskName[0] != '\0') {
			sprintf(strMessage,"Using user-defined image plane mask file <%s>",strMaskName);
			TRACE(LOG_HIGH_PRI,strMessage);
			pMask = lv_read_img_file(strMaskName,g_PixelResn);
			if (pMask==NULL){
				sprintf(strMessage,"Failed to read mask file <%s>.",strMaskName);
				LOG_ERR(strMessage);
				goto FAIL_EXIT;
			}
			if(pReadImg->pAxisSize[0] != pMask->pAxisSize[0] || pReadImg->pAxisSize[1] != pMask->pAxisSize[1]) {
				sprintf(strMessage,"ERROR: data image (%d,%d) and mask image (%d,%d) are not the same size",
						(int)pReadImg->pAxisSize[0],(int)pReadImg->pAxisSize[1],
						(int)pMask->pAxisSize[0],(int)pMask->pAxisSize[1]);
				LOG_ERR(strMessage);
				goto FAIL_EXIT;
			}
		}
	}

	/* make structures for the ray traced and map projected images */
	pProjImage = lv_create_image_struct(g_iSzImgx,g_iSzImgy, IMGTYPE, g_PixelResn);
	pMappedImage = lv_create_image_struct(g_iSzImgx,g_iSzImgy, IMGTYPE, g_PixelResn);

	/* copy the external mask into the mapped image struct mask */
	if (pMask!=NULL) {
		int iImgSize=0;

		iImgSize = img_CalcImgSize(pMask);
		pMappedImage->pMask = malloc(sizeof(bool)*iImgSize);
		if (pMappedImage->pMask==NULL) {
			LOG_ERR("No malloc for external image mask");
			goto FAIL_EXIT;
		}
		pMappedImage->bExternalMask=TRUE;
		for (i=0; i< iImgSize; i++) {
			pMappedImage->pMask[i] = (bool) pMask->pImage[i];
		}
	}

	/* read the PSF if it has been specified */
	if (strcmp(strPsfFilename,"")!= 0) {
		sprintf(strMessage,"Reading PSF file <%s>",strPsfFilename);
		TRACE(LOG_HIGH_PRI,strMessage);
		pPsfImg = lv_read_img_file(strPsfFilename,g_PixelResn);
		if (pPsfImg == NULL) {
			sprintf(strMessage,"Failed to read FITS file %s.",strPsfFilename);
			LOG_ERR(strMessage);
			fprintf(stderr,"%s\n",strMessage);
			goto EXIT;
		}

		/* normalise the PSF */
		if (iNormaliseMax !=0) {
			temp = img_MaxArrayVal(pPsfImg->pImage,img_CalcImgSize(pPsfImg));
			sprintf(strMessage,"Max value in PSF array was %g. Normalising so that max is 1...",temp);
		}
		else {
			temp = img_CalcArrayTotal(pPsfImg->pImage,img_CalcImgSize(pPsfImg));
			sprintf(strMessage,"Total of PSF array was %g. Normalising so that total is 1...",temp);
		}
		if ( fabs(temp-1.0) > 1e-5) {
			img_VecMultByScalar(pPsfImg->pImage,1/temp,img_CalcImgSize(pPsfImg));
			TRACE(LOG_HIGH_PRI,strMessage);
		}
	}

	if (fConstVar != 0) {
		sprintf(strMessage,"Using user-specified global variance value: %g",fConstVar);
		TRACE(LOG_HIGH_PRI,strMessage);
		g_imgImgVariance = fConstVar;
		g_imgCalcVarianceFlag = FALSE;
	}

	/* read the noise file if it has been specified */
	if (strcmp(strNoiseFileName,"")!= 0) {
		if (fConstVar != 0) {
			LOG_ERR("WARNING: user-defined constant variance AND variance file have been specified. Using variance file");
		}
		else {
			sprintf(strMessage,"Reading variance file <%s>",strNoiseFileName);
			TRACE(LOG_HIGH_PRI,strMessage);
		}
		pNoiseImg = lv_read_img_file(strNoiseFileName,g_PixelResn);
		if (pNoiseImg == NULL) {
			sprintf(strMessage,"Failed to read FITS file %s.",strNoiseFileName);
			LOG_ERR(strMessage);
			fprintf(stderr,"%s\n",strMessage);
			goto EXIT;
		}
		if (pNoiseImg->pAxisSize[0] != pReadImg->pAxisSize[0] || pNoiseImg->pAxisSize[1] != pReadImg->pAxisSize[1]) {
			sprintf(strMessage,"Data and noise images are not the same size. Exiting\n");
			LOG_ERR(strMessage);
			fprintf(stderr,"%s\n",strMessage);
			goto EXIT;
		}
		for (i=0; i< pNoiseImg->pAxisSize[0]*pNoiseImg->pAxisSize[1]; i++) {
			if (pNoiseImg->pImage[i] <= 0) {
				LOG_ERR("ERROR: Variance file has zero and/or negative values.");
				goto FAIL_EXIT;
			}
		}
		g_imgCalcVarianceFlag = TRUE;
	}
	if (strNoiseFileName[0] == '\0' && fConstVar==0) {
		LOG_ERR("WARNING: no user-specified variance file or constant variance value. Using image properties to calculate variance. This is dodgey.");
	}

	/* make a source or read one in */
	if (strcmp(strSourceFileName,"") == 0) {
		LOG_ERR("No source supplied");
		goto FAIL_EXIT;
	}
	else {
		sprintf(strMessage,"Reading source file <%s>. Including offsets x: %g, y: %g",strSourceFileName,x_axisoffset,y_axisoffset);
		TRACE(LOG_HIGH_PRI,strMessage);
		pSourceOriginal = lv_read_img_file(strSourceFileName,(fPixelResnRatio != 0.0 ? g_PixelResn/fPixelResnRatio : g_PixelResn));
		if (pSourceOriginal == NULL) {
			sprintf(strMessage,"Failed to read FITS file %s.",strSourceFileName);
			LOG_ERR(strMessage);
			fprintf(stderr,"%s\n",strMessage);
			goto EXIT;
		}
		pSourceOriginal->pAxisOffset[0] = x_axisoffset;
		pSourceOriginal->pAxisOffset[1] = y_axisoffset;
	}

/**** test convolution function
pTemp = lv_duplicate_image(pPsfImg);
img_ReverseTransposeImg(pPsfImg,pTemp);
newImg = img_ConvolveImgsWithFFT(pSourceOriginal,pTemp,NULL);
if (newImg == NULL) {
	fprintf(stderr,"newImg is NULL.\n");
}
else {
	lv_write_image_to_file(newImg,"test_convol_out.fits",TRUE);
}
goto EXIT;
*************/

	/* set a default value for the src plane if it hasn't
		been specified */
	if (fSrcDefaultVal == 0) {
		fSrcDefaultVal = img_CalcArrayTotal(pReadImg->pImage,img_CalcImgSize(pReadImg))/img_CalcImgSize(pSourceOriginal)/10.0;
		sprintf(strMessage,"Using calculated default src val of %g.",fSrcDefaultVal);
	}
	else {
		sprintf(strMessage,"Using user-specified default src val of %g.",fSrcDefaultVal);
	}
	TRACE(LOG_HIGH_PRI,strMessage);

	/* keep an original copy of the soure */
	pSource = lv_duplicate_image(pSourceOriginal);

	sprintf(strMessage,"Reading param file <%s>",strParamFile);
	TRACE(LOG_HIGH_PRI,strMessage);

	pLens = lm_ReadLensFile(strParamFile);
	if (pLens == NULL) {
		fprintf(stderr,"lm_ReadLensFile failed.\n");
		goto FAIL_EXIT;
	}

	sprintf(strMessage,"Image pixel resolution: %g. Source to image pixel resolution ratio: %g",g_PixelResn,fPixelResnRatio);
	TRACE(LOG_HIGH_PRI,strMessage)

    /* make an image using simple ray-tracing. no PSF */
	if (iRayTraceOnly ==1) {
		TRACE(LOG_HIGH_PRI,"Projecting image using Ray Tracing.");

		/* project the source using simple ray tracing */
		lv_projectImage(pLens, pSource, pProjImage);
	}

	/* create a mapping matrix for for LensMEM weights matrix */
	if (iRayTraceOnly ==0 && iMinimise == FALSE) {
	    /* project source through map matrix, no parameter searching */	
		real_t	tempchi=0;
		int	iNumImgPix=0;

	    /* allocate the mapping matrix array*/
	    pWeightMatrix = lv_allocMappingMatrix(pMappedImage->pAxisSize[0],pMappedImage->pAxisSize[1],pSource->pAxisSize[0],pSource->pAxisSize[1]);
	    if (pWeightMatrix == NULL) {
		    fprintf(stderr,"lv_allocMappingMatrix failed.\n");
		    goto EXIT;
	    }

		TRACE(LOG_HIGH_PRI,"Creating mapping matrix");
		lv_zeroMapMatrix(pWeightMatrix);
		iStatus = lv_createMappingMatrix(pLens, pWeightMatrix, pSource, pMappedImage);
		if (iStatus != 0) {
			fprintf(stderr,"lv_createMappingMatrix failed. Stopping...\n");
			goto FAIL_EXIT;
		}

		TRACE(LOG_HIGH_PRI,"Projecting source through mapping matrix");
		lv_projectSourceThruMapMatrix(pSource, pMappedImage, pWeightMatrix);
		if (pPsfImg != NULL) {
		    lv_image_t *newImg=NULL; 
			TRACE(LOG_HIGH_PRI,"Convolving mapped image with PSF.");
			newImg = img_ConvolveImgsWithFFT(pMappedImage,pPsfImg,NULL);
			tempchi = img_CalcChiSquared(pReadImg,newImg,0,pNoiseImg,pMappedImage->pMask);
			iNumImgPix = img_CountActiveImgPixels(pMappedImage);
			sprintf(strMessage,"chi-squ between read image and projected image is %g. Num active img pix: %d\n",tempchi,iNumImgPix);
			TRACE(LOG_HIGH_PRI,strMessage);
			lv_free_image_struct(pMappedImage);
			pMappedImage = newImg;
		}
	}

    if (iRayTraceOnly ==0 && iMakeSrcInvMagMap && (iMinimise==TRUE || iUseMinFinder==TRUE) ) {
        sprintf(strMessage,"ERROR: Cannot generate magnification map when parameter searching. Use a fixed set of params and no '-dofit' or '-useminfinder'.\n");
        TRACE(LOG_ERROR_PRI,strMessage);
        fprintf(stderr,strMessage);
        goto FAIL_EXIT;
    }

	if (iRayTraceOnly ==0 && iMinimise == TRUE && iUseMinFinder == FALSE) {
	    iStatus = lp_ParamSweep(pLens, pSource, pMappedImage, pReadImg,iMaxIterations, &fResChiSqu, &fResEntropy, &pBestSrc, &pBestImg, pPsfImg, pNoiseImg, iSearchCentre,fSearchCentreRange,fSearchCentreStep,fSrcDefaultVal,iMethod,iDumpImgs);
        if (iStatus !=0 ) {
	        sprintf(strMessage,"lp_ParamSweep returned %d. Exiting...",iStatus);
	        TRACE(LOG_ERROR_PRI,strMessage);
	        goto FAIL_EXIT;
        }
    }

	if (iRayTraceOnly ==0 && iUseMinFinder == TRUE) {
		iStatus = lp_MinFinder(pLens, pSource, pMappedImage, pReadImg,iMaxIterations, &fResChiSqu, &fResEntropy, &pBestSrc, &pBestImg, pPsfImg, pNoiseImg, iSearchCentre,fSrcDefaultVal);
        if (iStatus !=0 ) {
	        sprintf(strMessage,"lp_MinFinder returned %d. Exiting...",iStatus);
	        TRACE(LOG_ERROR_PRI,strMessage);
	        goto FAIL_EXIT;
        }
	}

	/* create img plane inverse mag, critline, src plane mag & caustic */
	if (iRayTraceOnly ==0 && iMakeSrcInvMagMap && iMinimise==FALSE && iUseMinFinder == FALSE) {
		TRACE(LOG_HIGH_PRI,"Creating and writing source plane inverse magnification map");
		pInvMag = img_CalcSrcPlaneMag(pWeightMatrix,pSourceOriginal->fPixelAngSize, g_PixelResn);
        lv_write_image_to_file(pInvMag,strOutSrcMagFilename,TRUE);
		for (i=0; i< img_CalcImgSize(pInvMag); i++) pInvMag->pImage[i] = pWeightMatrix->pSumSrc[i];
		lv_write_image_to_file(pInvMag,"src_sum.fits",TRUE);
		lv_free_image_struct(pInvMag);
		TRACE(LOG_HIGH_PRI,"Creating and writing image plane inverse mag map");
		pInvMag = img_CalcImgPlaneInvMag(pWeightMatrix,pSourceOriginal->fPixelAngSize, g_PixelResn);
		lv_write_image_to_file(pInvMag,strOutImgMagFilename,TRUE);
		/* find the critical line and the caustic */
		cl_CalcCritLinePoints(pLens,pInvMag,&cl_list_head);
		cl_CalcCausticPoints(pLens,&cl_list_head,&caust_list_head);
		/* print them to a file */
		TRACE(LOG_HIGH_PRI,"Dumping critical line and caustic points to text files.");
		cl_SavePixList(&cl_list_head,"critline.txt");
		cl_SavePixList(&caust_list_head,"caustic.txt");
		/* free */
		cl_FreePixList(cl_list_head.pNext);
		cl_FreePixList(caust_list_head.pNext);
		/* write a boolean map of the multiply imaged pixels */
		for(i=0; i<pInvMag->pAxisSize[0]*pInvMag->pAxisSize[1]; i++) {
			pInvMag->pImage[i] = (real_t) pWeightMatrix->pMultImgPix[i];
		}
		lv_write_image_to_file(pInvMag,"mult_img_pix.fits",TRUE);
		lv_free_image_struct(pInvMag); pInvMag=NULL;
		goto EXIT;
	}

	/* write resulting image to file */
	sprintf(strMessage,"Writing model images to file. Chisqu: %g, ent: %g",fResChiSqu,fResEntropy);
	TRACE(LOG_HIGH_PRI,strMessage);

	if (iRayTraceOnly ==0) {
	    /* with no minimisation, result is a projected image */
	    if (iRayTraceOnly ==0 && iMinimise == FALSE) {
	        //if (pMappedImage != NULL) lv_write_image_to_file(pMappedImage,strOutputMappedFileName, TRUE);
            if (pMappedImage != NULL) lv_write_image_to_file(pMappedImage,dynamicOutputModelImage, TRUE);
	    }
		//if (pBestImg != NULL) lv_write_image_to_file(pBestImg,strOutputMappedFileName, TRUE);
        if (pBestImg != NULL) {
            lv_write_image_to_file(pBestImg,dynamicOutputModelImage, TRUE);
            lv_image_t *residual = lv_residual_img_file(pReadImg, pBestImg);
            lv_write_image_to_file(residual, dynamicOutputResidulImage, TRUE);
            lv_free_image_struct(residual);
        }
        //if (pBestSrc != NULL) lv_write_image_to_file(pBestSrc,strOutputSourceFile,TRUE);
        if (pBestSrc != NULL) lv_write_image_to_file(pBestSrc,dynamicOutputSourceImage,TRUE);
        imgDumpMask(pMappedImage);
	}
	else {
		lv_write_image_to_file(pProjImage, strOutputImageFile, TRUE);
	}

FAIL_EXIT:

	TRACE(LOG_HIGH_PRI,"Freeing structures");
	if (pLens != NULL) lm_FreeLensModel(pLens);
	if (pTemp != NULL) lv_free_image_struct(pTemp);
	if (pSource != NULL) lv_free_image_struct(pSource);
	if (pBestSrc != NULL) lv_free_image_struct(pBestSrc);
	if (pBestImg != NULL) lv_free_image_struct(pBestImg);
	if (pSourceOriginal != NULL) lv_free_image_struct(pSourceOriginal);
	if (pPsfImg != NULL) lv_free_image_struct(pPsfImg);
	if (pReadImg != NULL) lv_free_image_struct(pReadImg);
	if (pNoiseImg != NULL) lv_free_image_struct(pNoiseImg);
	if (pMask != NULL) lv_free_image_struct(pMask);
	if (pMappedImage != NULL) lv_free_image_struct(pMappedImage);
	if (pProjImage != NULL) lv_free_image_struct(pProjImage);
	if (pInvMag != NULL) lv_free_image_struct(pInvMag);
	if (pWeightMatrix != NULL) lv_freeMappingMatrix(pWeightMatrix);

EXIT:

	TRACE_OUT;
	return 0;
}

/******************************
 *****************************/
void printWeight(int nrow, int ncol, float *pWeight) {
	int i,j;

	for(j=0; j < nrow; j++) {
		for (i=0; i< ncol; i++) {
			printf("%f, ",*pWeight);
			pWeight++;
		}
		printf("\n");
	}
}


/***************************
 *************************/
void printUsage(int argc, char * const argv[]) {
  fprintf(stderr,"Usage:\n");
  fprintf(stderr,"-imgfile filename\n");
  fprintf(stderr,"-logfilepath filepath\n");
  fprintf(stderr,"-srcxoffset offset\n");
  fprintf(stderr,"-srcyoffset offset\n");
  fprintf(stderr,"-tracelevel level : logging level. 3=essential messages, 2=debugging messages, 1=huge amount of debugging messages\n");
  fprintf(stderr,"-raytraceonly : flag to just ray-trace the source, make an image and exit\n");
  fprintf(stderr,"-sourcefile filename\n");
  fprintf(stderr,"-psffile filename\n");
  fprintf(stderr,"-datafile filename\n");
  fprintf(stderr,"-dofit : flag to turn on parameter fitting\n");
  fprintf(stderr,"-nice num : equivalent to Unix nice\n");
  fprintf(stderr,"-pixelratio ratio : =(angular size of pixels in image plane)/(angular size of pixels in source plane)\n");
  fprintf(stderr,"-maxiter num : maximum number of iterations in MEM loop \n");
  fprintf(stderr,"-pixelres num : angular size of pixels in image (arcsec)\n");
  fprintf(stderr,"-noisefile filename\n");
  fprintf(stderr,"-makemag : flag to make image and source plane magnification maps and trace the critcal line and caustic (experimental)\n");
  fprintf(stderr,"-dumpimg : flag to have (some) images generated by parameter search written out\n");
  fprintf(stderr,"-paramfile filename\n");
  fprintf(stderr,"-fixlambda num\n");
  fprintf(stderr,"-srcdefaultval num\n");
  fprintf(stderr,"-mask filename\n");
  fprintf(stderr,"-useconjgrad : (obsolete) flag to use conjugate gradient minimisation\n");
  fprintf(stderr,"-normalisepsfmax : flag to have the PSF file normalised so that its maximum value is 1.0\n");
  fprintf(stderr,"-debugimgs : produce huge number of images used for debugging MEM loop\n");
  fprintf(stderr,"-targetchisqu num\n");
  fprintf(stderr,"-usemultimgpix : (experimental) \n");
  fprintf(stderr,"-constvariance num : specify a fixed variance for all image pixels\n");
  fprintf(stderr,"-searchcentre : (experimental) flag to enable the lens centre to be free parameters\n");
  fprintf(stderr,"-searchcentrerange num : range of pixels for lens centre as free param\n");
  fprintf(stderr,"-searchcentrestep num : step size in pixels for lens centre as free param\n");
  fprintf(stderr,"-useminfinder : flag to use an amoeba simplex to search for best fit params\n");

  exit(0);
}
