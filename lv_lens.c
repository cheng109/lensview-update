/*******************
 	functions for creating lens models and
	calculating deflections etc of lens models
	for the "lensview" software
Copyright (C) 2006. Randall Wayth.
$Id: lv_lens.c,v 1.32 2008/10/31 20:43:10 rwayth Exp rwayth $
*******************/
#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<string.h>
#include	"common.h"
#include	"log.h"
#include	"lv_common.h"
#include	"lv_lens.h"
#include	"fastell.h"
#include	"gsl/gsl_integration.h"

#define		BUFF_SIZE	256
/* NOTE: if any lens model has > 5 params, then MAX_LENSCOMP_PARAMS
 * will need to be changed in lv_common.h */
#define	LM_NPARAMS_EXTSH	2
#define	LM_NPARAMS_NFW		4
#define	LM_NPARAMS_PTMASS	1
#define	LM_NPARAMS_MASSSHEET	1
#define	LM_NPARAMS_SPEMD	7
#define	LM_NPARAMS_SIE		3
#define	LM_NPARAMS_PIEP		4
#define	LM_NPARAMS_EXPDISC	4
#define	LM_NPARAMS_SIS		1
#define	LM_NPARAMS_SERSIC	5
#define	LM_NPARAMS_FERRERS	4
#define	LM_NPARAMS_USERDEF	2

#define LM_NPARAMS_SPEMD_NEW    7 

#define	GOLD_RATIO	1.618034

#define	LM_CACHE_MISS		-42
typedef	struct	_lmDeflCache {
	short			iType;
	real_t			param[3];		/* axratio, scale, r_power */
	real_t			pix_size;
	lv_axissize_t	dimension[2];
	real_t			*pValsX;
	real_t			*pValsY;
	real_t			*pDeflX;
	real_t			*pDeflY;
	struct  _lmDeflCache *pNext;
}	lmDeflCache;

static	double	lm_arctanh(double x);
static	double	lm_arccosh(double x);
static	double	lm_nfw_mass(double x);
/*static	real_t	lm_expdiscintx(real_t u);*/
/*static	real_t	lm_expdiscinty(real_t u);*/
static  double  lm_expdiscintx_gsl(double x, void *params);
static  double  lm_expdiscinty_gsl(double x, void *params);
static	real_t	lm_hernqf(real_t x);
static	real_t	lm_CalcDeVaucCircDefl(real_t r, real_t Re);
/*static	real_t	lm_SersicDeflX(real_t x);
static	real_t	lm_SersicDeflY(real_t x);*/
static	double	lm_SersicDeflX_gsl(double x, void *params);
static	double	lm_SersicDeflY_gsl(double x, void *params);
static	int 	lm_deflCacheLookup(short iType, lmDeflCache *pCache, real_t x, real_t y, real_t fAxratio, real_t fScale, real_t fM, real_t *pDeflx, real_t *pDefly);
static  real_t   lm_calc_sersic_b(real_t m);
static int	lm_CreateLMComp_SPEMD(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t fMassScaleFrom, 
		real_t fMassScaleTo, real_t fMassScaleInc, real_t fEllipFrom, real_t fEllipTo, real_t fEllipInc,
		real_t fAngleFrom, real_t fAngleTo, real_t fAngleInc, real_t fGammaFrom, real_t fGammaTo,
		real_t fGammaInc, real_t fCoreFrom, real_t fCoreTo, real_t fCoreInc,
		real_t fCenterXFromm, real_t fCenterXTo, real_t fCenterXInc,
		real_t fCenterYFrom, real_t CenterYTo, real_t CenterYInc);
static int	lm_CreateLMComp_Ferrers(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t critfrom, real_t critto, real_t critinc, real_t AFrom, real_t ATo, real_t AInc, real_t BFrom, real_t BTo, real_t BInc, real_t Angfrom, real_t Angto, real_t Anginc);
static int	lm_CreateLMComp_Userdef(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t critfrom, real_t critto, real_t critinc, real_t AFrom, real_t ATo, real_t AInc, lv_image_t *pDefX, lv_image_t *pDefY);
static	double	lm_ferrersq10(real_t a, real_t b, double lambda);
static	double	lm_ferrersq01(real_t a, real_t b, double lambda);
static	double	lm_ferrersq03(real_t a_squ, real_t b_squ, double l);
static	double	lm_ferrersq30(real_t a_squ, real_t b_squ, double l);
static	double	lm_ferrersq12(real_t a_squ, real_t b_squ, double l);
static	double	lm_ferrersq21(real_t a_squ, real_t b_squ, double l);
static	double	lm_ferrersq02(real_t a_squ, real_t b_squ, double l);
static	double	lm_ferrersq20(real_t a_squ, real_t b_squ, double l);
static	double	lm_ferrersq11(real_t a_squ, real_t b_squ, double lambda);
static	int	lm_CalcFerrersDefl(real_t x, real_t y, real_t a, real_t b, real_t *fResX, real_t *fResY);

/* external global variables */
extern	lv_axissize_t	g_iSzImgx, g_iSzImgy;
float  g_PixelResn= HST_WF_RESOLUTION/2.0;

/* private global variables */
static	real_t	g_tempx,g_tempy,g_scale,g_axratio;
static	real_t	g_sersic_b;
static	int	g_iCompNum=0;

/***************************
Function:	lm_CalcDeflection
Description:	Calculate the deflection of a ray in polar coords normalised to Critical radius units.
				Deflection is relative to the centre of the lens. Lens offsets must be calculatd
				before and after this.
Arguments:	pLens:	Lens which is doing the deflecting. Assumes initialised etc.
			fX, fY: coords (x,y) of ray in image plane in angular units
			fDeltaY, pDeltaY: deflection vector (x,y)
			pMagnification: mag caused by deflection. Not implemented
Returns:	0: success
			other: failure
****************************/
int	lm_CalcDeflection(lv_lensmodel_t *pLens, real_t fX, real_t fY, real_t *pDeltaX, real_t *pDeltaY, real_t *pMagnification) {
	int	i,iStatus=0;
	real_t	delxcomp =0, delycomp=0;

	for (i=0; i< pLens->iNumComponents; i++) {
		delxcomp =0; delycomp=0;
		iStatus = lm_CalcDeflComponent(pLens->pComponent + i, fX, fY, &delxcomp, &delycomp);
/*		printf("deflection at (%g,%g) is (%g,%g)\n",fX,fY,delxcomp,delycomp);  */
		if (iStatus != 0) {
			goto EXIT;
		}
		*pDeltaX += delxcomp;
		*pDeltaY += delycomp;
	}

EXIT:
	/* calculate the magnification */
	/******* NOT IMPLEMENTED *********/
	*pMagnification = 1.0;
	return iStatus;
}


/***************************
Function:		lm_CalcDeflComponent
Description:
Arguments:		fX, fY: coordinate in image plane in angle (arcsec) units.
Returns:
****************************/
int lm_CalcDeflComponent(lv_lenscomp *pLensComp, real_t fX, real_t fY, real_t *pDeltaX, real_t *pDeltaY){

	static	lmDeflCache	deVaucCache = {-1,{-1,-1,-1},0,{0,0},NULL,NULL,NULL,NULL,NULL};
	static	lmDeflCache	expDiscCache = {-1,{-1,-1,-1},0,{0,0},NULL,NULL,NULL,NULL,NULL};
	int		iStatus =0;

	/* add the offset for each lens component */
	fX -= g_PixelResn*pLensComp->fXoffset;
	fY -= g_PixelResn*pLensComp->fYoffset;

	switch(pLensComp->iType) {
		case	LM_ISO_SPHERE:	 /* calculate for Isothermal sphere */
			{
				real_t	fDenom=0;

				fDenom = sqrt(fX*fX + fY*fY);
				if (fDenom == 0.0) {
					iStatus =LM_IGNORE_PROJ;
					*pDeltaX = pLensComp->fParameter[0];
					*pDeltaY = pLensComp->fParameter[0];
				}
				else {
					*pDeltaX = pLensComp->fParameter[0]*(fX/fDenom);
					*pDeltaY = pLensComp->fParameter[0]*(fY/fDenom);
				}
			}
			break;

		case	LM_PIEP:
			/* parameters: mass, ellipticity, angle, core */
			{
				real_t	fCritRad,fDenom,fCoreRad,fEllip,fCosTheta_g,fSinTheta_g;
				real_t	deltax1,deltay1,x1,y1;

				/* pre-calculate constants */
				fCosTheta_g = cos(pLensComp->fParameter[2]*M_PI/180);
				fSinTheta_g = sin(pLensComp->fParameter[2]*M_PI/180);

				fCritRad = pLensComp->fParameter[0];
				fCoreRad = pLensComp->fParameter[3];
				fEllip = pLensComp->fParameter[1];
                if (fEllip >= 1.0 || fCoreRad < 0 || fCritRad < 0) {
                    iStatus = LM_BADPROJ;
                    break;
                }

				/* rotate the reference frame for ellipse with semi-major axis of potential in X dir */
				x1 = fX*fCosTheta_g + fY*fSinTheta_g;
				y1 = -fX*fSinTheta_g + fY*fCosTheta_g;

				fDenom = sqrt(fCoreRad*fCoreRad + (1.0-fEllip)*x1*x1 + (1.0+fEllip)*y1*y1);

				if (fDenom != 0.0) {
					deltax1 = fCritRad*x1*(1.0-fEllip)/fDenom;
					deltay1 = fCritRad*y1*(1.0+fEllip)/fDenom;

					/* rotate back again */
					*pDeltaX = deltax1*fCosTheta_g - deltay1*fSinTheta_g;
					*pDeltaY = deltay1*fCosTheta_g + deltax1*fSinTheta_g;
				}
				else {
					/* for the dead centre of the lens plane, just deflect by the critical radius */
					iStatus =LM_IGNORE_PROJ;
					*pDeltaX = fCritRad*fCosTheta_g;
					*pDeltaY = fCritRad*fSinTheta_g;
				}
				break;
			}
		case	LM_PTMASS:
			{
				real_t	fDenom =0, fMult=0;
				/* NOTE this is r squared. */
				fDenom = (fX*fX + fY*fY);

				if (fDenom == 0.0) {
					/* for the dead centre of the lens plane, just deflect by a lot */
					*pDeltaX = 1000;
					*pDeltaY = 1000;
					iStatus =LM_IGNORE_PROJ;
				}
				else {

					fMult = (pLensComp->fParameter[0]*pLensComp->fParameter[0])/fDenom;

					/* calculate the deflection in x,y co-ords */
					*pDeltaX = fX*fMult;
					*pDeltaY = fY*fMult;
				}
			}
			break;

		case	LM_SPEMD:
			{
				/* parameters are: kappa, ellipticity,angle, gamma, corerad, centerX, centerY */

			  double	fTempKappa = 0.0,fTempCoreSqu=0.0, fTempAxratio=1.0, fTempDefl[2], fTempGamma = 0.0, fTempCenterX=0, fTempCenterY=0;
				double	x1, y1;
				real_t	fCosTheta, fSinTheta;

				/* note: need the sqrt(fTempAxratio) to normalise the area pi*a*b equal to one
				 * for all ellipticities */
				fTempAxratio = 1.0 - pLensComp->fParameter[1];
				fTempCoreSqu = pLensComp->fParameter[4]*pLensComp->fParameter[4];
				fTempGamma = pLensComp->fParameter[3];
                if (fTempAxratio > 1.0 || fTempAxratio<=0 ) {
                    iStatus = LM_BADPROJ;
                    break;
                }

				/*  if this is the centre of the lens plane, then the deflection is either
				 * 	zero for gamma < 0.5, and very large for gamma > 0.5. For gamma == 0.5
				 * (isothermal) then the deflection is the critical radius which is good
				 * enough in most cases since this will miss the source by miles */
				if (fX == 0 && fY == 0 && fTempGamma >= 0.5) {
					*pDeltaX = *pDeltaY = pLensComp->fParameter[0];
					iStatus =LM_IGNORE_PROJ;
					break;
				}

				/* pre-calculate constants */
				fCosTheta = cos(pLensComp->fParameter[2]*M_PI/180);
				fSinTheta = sin(pLensComp->fParameter[2]*M_PI/180);

				/* rotate the reference frame for ellipse with semi-major axis of potential in X dir */
				x1 = (fX*fCosTheta + fY*fSinTheta);
				y1 = (-fX*fSinTheta + fY*fCosTheta);

/* 20/6/01. Correction for non-isothermal added. See notebook
*				fTempKappa = 0.5*pLensComp->fParameter[0]/sqrt(fTempAxratio);
*				fTempKappa = 0.5*pLensComp->fParameter[0]/(pow(fTempAxratio,fTempGamma)*pow(pLensComp->fParameter[0],2.0*fTempGamma-1))*(2.0*fTempGamma);
*				fTempKappa = 0.5*pLensComp->fParameter[0]/(pow(fTempAxratio,fTempGamma)*(2.0-2.0*fTempGamma));
*				fTempKappa = 0.5*pLensComp->fParameter[0]/(pow(fTempAxratio,fTempGamma))*(2.0*fTempGamma);
				E=(2.0*(1-gam)/ellip)^(gam)
*/

				fTempKappa = 0.5*pLensComp->fParameter[0]*pow((2.0-2.0*fTempGamma)/fTempAxratio,fTempGamma);
				/* call the FORTRAN function thing */
				fastelldefl_(&x1,&y1,&fTempKappa,&fTempGamma,&fTempAxratio,&fTempCoreSqu,fTempDefl);

				/* rotate back again */
				*pDeltaX = fTempDefl[0]*fCosTheta - fTempDefl[1]*fSinTheta;
				*pDeltaY = fTempDefl[1]*fCosTheta + fTempDefl[0]*fSinTheta;
			}
			break;

		case	LM_NFW:
			/* parameters:  mass scale, scale len, ellipticity, orientation_angle */
			{
				real_t	fEllip,fCosTheta,fSinTheta,x1,y1,fPhi,fAngRadius,fTempResult,fCosPhi,fSinPhi,fScale;

				/* pre-calculate constants */
				fCosTheta = cos(pLensComp->fParameter[3]*M_PI/180);
				fSinTheta = sin(pLensComp->fParameter[3]*M_PI/180);

				fEllip = pLensComp->fParameter[2];
				fScale = pLensComp->fParameter[1];

                if (fEllip >= 1.0 || fEllip < 0 || fScale < 0) {
                    iStatus = LM_BADPROJ;
                    break;
                }

				/* create elliptical co-ords still in angle units from rotated frame sky coords */
				x1 = sqrt(1.0 - fEllip)*(fX*fCosTheta + fY*fSinTheta);
				y1 = sqrt(1.0 + fEllip)*(-fX*fSinTheta + fCosTheta*fY);
				fPhi = atan2(y1,x1);

				/* angular radius is in dimensionless units */
				fAngRadius = sqrt(x1*x1 + y1*y1)/fScale;

				if (fAngRadius > 0.0) {
					real_t	deflx,defly;

					fCosPhi = cos(fPhi);
					fSinPhi = sin(fPhi);
					/* calculate mass. Note that one multiple of fAngRadius has been removed from denominator in next
					 * line so that we don't need to mult again in the following two */
					fTempResult = pLensComp->fParameter[0]*fScale*lm_nfw_mass(fAngRadius)/(fAngRadius);
					deflx = sqrt(1.-fEllip)*fTempResult*fCosPhi;
					defly = sqrt(1.+fEllip)*fTempResult*fSinPhi;
					*pDeltaX = (deflx*fCosTheta - defly*fSinTheta);
					*pDeltaY = (deflx*fSinTheta + defly*fCosTheta);
				}
				else {
/*
					iStatus =LM_IGNORE_PROJ;
*/
					*pDeltaX = 0.0;
					*pDeltaY = 0.0;
				}
			}
			break;

		case	LM_SIE:
			/* parameters are: critrad, axratio, angle */
			{
				real_t	phi,root1mq,fq,fac,fCore=0,fCosTheta,fSinTheta,x1,y1,deltax1,deltay1;

				if (fX == 0 && fY == 0) {
					*pDeltaX = *pDeltaY = pLensComp->fParameter[0];
					iStatus =LM_IGNORE_PROJ;
					break;
				}

				/* pre-calculate constants */
				fCosTheta = cos(pLensComp->fParameter[2]*M_PI/180);
				fSinTheta = sin(pLensComp->fParameter[2]*M_PI/180);

				fq = pLensComp->fParameter[1];
                if (fq > 1.0) {
                    iStatus = LM_BADPROJ;
                    break;
                }
                if (fq==1.0) fq = 0.999;

				/* rotate reference frame to x-axis */
				x1 = fX*fCosTheta + fY*fSinTheta;
				y1 = -fX*fSinTheta + fY*fCosTheta;

				root1mq = sqrt(1.0-fq*fq);
				phi = sqrt(fq*fq*(fCore*fCore + x1*x1) + y1*y1);
				/* use sqrt(fq) here not fq. This presevers scale of Einstein Radius */
				fac = pLensComp->fParameter[0]*sqrt(fq)/root1mq;

				deltax1 = fac*atan(root1mq*x1/(phi + fCore));
				deltay1 = fac*lm_arctanh(root1mq*y1/(phi+ fCore*fq*fq));

				/* rotate back again */
				*pDeltaX = deltax1*fCosTheta - deltay1*fSinTheta;
				*pDeltaY = deltay1*fCosTheta + deltax1*fSinTheta;

                /* printf("x1,y1: %g,%g. root1mq: %g, phi: %g, fac: %g\n",x1,y1,root1mq,phi,fac); */

			}
			break;

		case	LM_EXTSHEAR:
			/* parameters: shear, shear angle */
			{
				real_t	x1,y1,fCosTheta,fSinTheta;

				fCosTheta = cos(pLensComp->fParameter[1]*M_PI/180);
				fSinTheta = sin(pLensComp->fParameter[1]*M_PI/180);

				x1 = fX*fCosTheta + fY*fSinTheta;
				y1 = -fX*fSinTheta + fY*fCosTheta;

				/* rotate back again */
				*pDeltaX = (x1*fCosTheta + y1*fSinTheta)* pLensComp->fParameter[0];
				*pDeltaY = (-y1*fCosTheta + x1*fSinTheta)* pLensComp->fParameter[0];

/*
				*pDeltaX = 2.0*fShear*fCosTheta*(fX*fCosTheta + fY*fSinTheta) - fShear*fX;
				*pDeltaY = 2.0*fShear*fSinTheta*(fX*fCosTheta + fY*fSinTheta) - fShear*fY;
*/
			}
			break;

		case	LM_MASSSHEET:
			{
				*pDeltaX = 2*pLensComp->fParameter[0]*fX;
				*pDeltaY = 2*pLensComp->fParameter[0]*fY;
			}
			break;

		case	LM_EXPDISC:
			/* parameters: kappa,inclination axis ratio (0-1),angle,scale len */
			{
				real_t	fCosTheta,fSinTheta,x1,y1,deflx,defly;
                if (pLensComp->fParameter[1] > 1.0 || pLensComp->fParameter[1] <= 0) {
                    iStatus = LM_BADPROJ;
                    break;
                }

				/* rotate so that major axis of disc in x direction */
				fCosTheta = cos(pLensComp->fParameter[2]*M_PI/180);
				fSinTheta = sin(pLensComp->fParameter[2]*M_PI/180);
				x1 = fX*fCosTheta + fY*fSinTheta;
				y1 = -fX*fSinTheta + fY*fCosTheta;

				iStatus = lm_deflCacheLookup(LM_EXPDISC,&expDiscCache,x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[3],0, &deflx, &defly);
				if (iStatus == LM_CACHE_MISS) {
					sprintf(strMessage,"WARNING: Cache miss for EXP DISC");
					deflx = lm_expdiscdeflx(x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[3]);
					defly = lm_expdiscdefly(x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[3]);
					iStatus = 0;
				}

				/* rotate back again */
				*pDeltaX = (deflx*fCosTheta - defly*fSinTheta)* pLensComp->fParameter[0];
				*pDeltaY = (defly*fCosTheta + deflx*fSinTheta)* pLensComp->fParameter[0];
			}
			break;
		case	LM_SERSIC:
			/* parameters: kappa,axis ratio (0-1),angle,scale len, M parameter */
			{
				real_t	fCosTheta,fSinTheta,x1,y1,deflx,defly;

                if (pLensComp->fParameter[1] > 1.0 || pLensComp->fParameter[1] <= 0) {
                    iStatus = LM_BADPROJ;
                    break;
                }
				/* rotate so that major axis of mass in x direction */
				fCosTheta = cos(pLensComp->fParameter[2]*M_PI/180);
				fSinTheta = sin(pLensComp->fParameter[2]*M_PI/180);
				x1 = fX*fCosTheta + fY*fSinTheta;
				y1 = -fX*fSinTheta + fY*fCosTheta;

				iStatus = lm_deflCacheLookup(LM_SERSIC,&deVaucCache,x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[3],pLensComp->fParameter[4], &deflx, &defly);
				if (iStatus == LM_CACHE_MISS) {
					iStatus = lm_CalcSersicDefl(x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[3],pLensComp->fParameter[4], &deflx, &defly);
					if (iStatus != 0) {
						goto EXIT;
					}
				}

				/* rotate back again */
				*pDeltaX = (deflx*fCosTheta - defly*fSinTheta)*pLensComp->fParameter[0];
				*pDeltaY = (defly*fCosTheta + deflx*fSinTheta)*pLensComp->fParameter[0];

			}
			break;

		case	LM_FERRERS:
			/* parameters: kappa,axis a, axis b, angle */
			{
				real_t	fCosTheta,fSinTheta,x1,y1,deflx,defly;

				/* rotate so that major axis of mass in x direction */
				fCosTheta = cos(pLensComp->fParameter[3]*M_PI/180);
				fSinTheta = sin(pLensComp->fParameter[3]*M_PI/180);
				x1 = fX*fCosTheta + fY*fSinTheta;
				y1 = -fX*fSinTheta + fY*fCosTheta;

				iStatus = lm_CalcFerrersDefl(x1,y1,pLensComp->fParameter[1],pLensComp->fParameter[2], &deflx, &defly);

				/* rotate back again */
				*pDeltaX = (deflx*fCosTheta - defly*fSinTheta)*pLensComp->fParameter[0];
				*pDeltaY = (defly*fCosTheta + deflx*fSinTheta)*pLensComp->fParameter[0];

			}
			break;
				
		case	LM_USERDEF:
			/* parameters: kappa,angle */
			{
				real_t	fCosTheta,fSinTheta,x1,y1,deflx,defly;
				real_t	t,u;
				lv_axissize_t	sizx,sizy,xu,xl,yu,yl;

				/* rotate so that major axis of mass in x direction */
				fCosTheta = cos(pLensComp->fParameter[1]*M_PI/180);
				fSinTheta = sin(pLensComp->fParameter[1]*M_PI/180);
				sizx = pLensComp->pDeflX->pAxisSize[0];
				sizy = pLensComp->pDeflY->pAxisSize[1];
				x1 = (fX*fCosTheta + fY*fSinTheta)/pLensComp->pDeflX->fPixelAngSize + sizx/2;;
				y1 = (-fX*fSinTheta + fY*fCosTheta)/pLensComp->pDeflY->fPixelAngSize + sizy/2;;

				/* we want to interpolate. Get pixels which are closest to desired point */
				xu = ceil(x1);			/* upper bound on pixel value */
				xl = floor(x1);			/* lower bound */
				if (xu >= sizx) xu = sizx-1;
				if (xl < 0) xl = 0;

				yu = ceil(y1);
				yl = floor(y1);
				if (yu >= sizy) yu = sizy -1;
				if (yl < 0) yl =0;

				if (xl == xu && xu < sizx-1) {
					xu += 1;
				}
				else if (xl==xu) {
					xl-=1;
				}

				if (yl == yu && yu < sizy-1){
					yu += 1;
				}
				else if(yl == yu){
					yl-=1;
				}

				/* interpolate between grid points */
				t = (x1-xl)/(real_t)(xu-xl);
				u = (y1-yl)/(real_t)(yu-yl);
				if (isnan(t)|| isnan(u)) fprintf(stderr,"nan for x1,y1 (%g,%g), x: (%d,%d), y: (%d,%d) t,u: (%g,%g)\n",
								x1,y1,(int)xu,(int)xl,(int)yu,(int)yl,t,u);
				deflx = (1.0-t)*(1.0-u)*pLensComp->pDeflX->pImage[yl*sizx+xl]
						+ t*(1.0-u)*pLensComp->pDeflX->pImage[yl*sizx+xu]
						+ t*u*pLensComp->pDeflX->pImage[yu*sizx+xu]
						+ (1.0-t)*u*pLensComp->pDeflX->pImage[yu*sizx+xl];
				if (isnan(deflx)){
						fprintf(stderr,"deflx nan for x1,y1 (%g,%g), x: (%d,%d), y: (%d,%d) t,u: (%g,%g)\n",
								x1,y1,(int)xu,(int)xl,(int)yu,(int)yl,t,u);
						fprintf(stderr,"deflxs: %g,%g,%g,%g\n",pLensComp->pDeflX->pImage[yl*sizx+xl],
										pLensComp->pDeflX->pImage[yl*sizx+xu],
										pLensComp->pDeflX->pImage[yu*sizx+xu],
										pLensComp->pDeflX->pImage[yu*sizx+xl]);
				}
				defly = (1.0-t)*(1.0-u)*pLensComp->pDeflY->pImage[yl*sizy+xl]
						+ t*(1.0-u)*pLensComp->pDeflY->pImage[yl*sizy+xu]
						+ t*u*pLensComp->pDeflY->pImage[yu*sizy+xu]
						+ (1.0-t)*u*pLensComp->pDeflY->pImage[yu*sizy+xl];
				if (isnan(defly)){
						fprintf(stderr,"defly nan for x1,y1 (%g,%g), x: (%d,%d), y: (%d,%d) t,u: (%g,%g)\n",
								x1,y1,(int)xu,(int)xl,(int)yu,(int)yl,t,u);
						fprintf(stderr,"deflys: %g,%g,%g,%g\n",
										pLensComp->pDeflY->pImage[yl*sizy+xl],
										pLensComp->pDeflY->pImage[yl*sizy+xu],
										pLensComp->pDeflY->pImage[yu*sizy+xu],
										pLensComp->pDeflY->pImage[yu*sizy+xl]);
				}
				/* printf("x,y: %g, %g. x1,y1: %g, %g. defx: %g, defy: %g\n",fX,fY,x1,y1,deflx,defly); */

				/* rotate back again */
				*pDeltaX = (deflx*fCosTheta - defly*fSinTheta)*pLensComp->fParameter[0];
				*pDeltaY = (defly*fCosTheta + deflx*fSinTheta)*pLensComp->fParameter[0];

			}
			break;
		default:
			fprintf(stderr,"lm_CalcDeflComponent: Unknown lens model type %d. Returning failure.\n",pLensComp->iType);
			iStatus = 1;
			break;
	}
    if (isnan(*pDeltaX) || isnan(*pDeltaY) || isinf(*pDeltaX) || isinf(*pDeltaY) ) {
        iStatus =LM_BADPROJ;
        //fprintf(stderr,"ERROR: nan detected in lm_CalcDeflComponent\n");
    }

EXIT:
	return iStatus;
}


/***************************
Function:	lm_ReadLensFile
Description:	Read the lens components from a file
Arguments:
Returns:
****************************/
lv_lensmodel_t *lm_ReadLensFile(char	*strFile) {
	FILE	*fp=NULL;
	lv_lensmodel_t	*pLens= NULL;
	int		iNumComp=0, nconv=0, iCompType=0,iStatus=0;
	char	buffer[BUFF_SIZE],*pConv=NULL,strCompName[100],*pStartParam=NULL;

	TRACE_IN(lm_ReadLensFile);
	if (strFile == NULL) {
		fprintf(stderr,"lm_ReadLensFile: file name was null\n");
		goto EXIT;
	}

	if ((fp=fopen(strFile,"r")) == NULL) {
		sprintf(strMessage,"Could not open file <%s>",strFile);
		LOG_ERR(strMessage);
		goto EXIT;
	}

	/* count the number of lens components in the file */
	while (!feof(fp)) {
		pConv = fgets(buffer,BUFF_SIZE,fp);
		if (pConv==NULL) {
			break;
		}
		/* lines starting with '#' are comments. Skip comments and blanks */
		if (buffer[0] != '#' && buffer[0] != '\n') {
			iNumComp++;
		}
	}

	sprintf(strMessage,"Found %d lens component%s.",iNumComp,iNumComp==1? "" : "s");
	TRACE(LOG_MED_PRI,strMessage);

	if (iNumComp == 0) {
		goto EXIT;
	}

	rewind(fp);

	/* init a lens model */
	pLens = lm_CreateLensModel(iNumComp);

	if (pLens == NULL) {
		TRACE(LOG_ERROR_PRI,"failed to create lens model")
		goto EXIT;
	}

	/* now read each lens component and initialise it */
	while (!feof(fp)) {
		pConv = fgets(buffer,BUFF_SIZE,fp);
		if (pConv==NULL) {
			break;
		}
		/* skip comment/blank lines */
		if (buffer[0] != '#' && buffer[0] != '\n') {

			pStartParam = strstr(buffer,LM_RHMARK);

			if (pStartParam == NULL) {
				fprintf(stderr,"unexpected end of line. Can't find LM_RHMARK in input line <%s>\n",buffer);
				continue;
			}

			/* replace closing marker with whitespace so the scanf below will stop at the end
			 * of the component type name */
			*pStartParam = ' ';

			/* read the start of the line which identifies the component type */
			nconv = sscanf(buffer,LM_NM_LENSCOMP LM_LHMARK "%s" ,strCompName);
			if (nconv < 1) {
				fprintf(stderr,"Failed to scanf line <%s> for lens comp name\n",buffer);
				continue;
			}

			/* decode the type and set up that component*/
			iCompType = lm_decodeLensCompType(strCompName);
			if (iCompType <= 0) {
				sprintf(strMessage,"unknown lens comp type <%s>",strCompName);
				TRACE(LOG_ERROR_PRI,strMessage);
				continue;
			}

/*
			printf("Found lens comp %s type %d\n",strCompName,iCompType);
*/

			/* now read the rest of the line depending on the model */
			switch(iCompType) {
				case	LM_PIEP:
					{
						double	xoff=0,yoff=0,critfrom=0,critto=0,critinc=0,ellfrom=0,ellto=0,ellinc=0;
						double	angfrom=0,angto=0,anginc=0,corefrom=0,coreto=0,coreinc=0;

						nconv = sscanf(pStartParam+1,"(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf)",
                                       &xoff,&yoff,&critfrom,&critto,&critinc,
                                       &ellfrom,&ellto,&ellinc,&angfrom,&angto,&anginc,&corefrom,&coreto,&coreinc);

						if (nconv < LM_NPARAMS_PIEP*3+2) {
							fprintf(stderr,"Parameter conversion failed for lens type %d. params were: <%s>\n",iCompType,pStartParam+1);
							iStatus = 1;
							break;
						}

						iStatus = lm_CreateLMComp_PIEP(pLens,xoff,yoff,critfrom,critto,critinc,ellfrom,ellto,ellinc,
							angfrom,angto,anginc,corefrom,coreto,coreinc);

					}
					break;
				case	LM_EXTSHEAR:
					{
						double	shfrom=0,shto=0,shinc=0,angfrom=0,angto=0,anginc=0;

						nconv = sscanf(pStartParam+1,"(%lf,%lf,%lf,%lf,%lf,%lf)",&shfrom,&shto,&shinc,&angfrom, &angto,&anginc);
						if (nconv < LM_NPARAMS_EXTSH*3) {
							fprintf(stderr,"Parameter conversion failed for lens type %d. params were: <%s>\n",iCompType,pStartParam+1);
							iStatus = 1;
							break;
						}

						iStatus = lm_CreateLMComp_ExtShear(pLens,shfrom,shto,shinc,angfrom,angto,anginc);
					}
					break;
				case	LM_MASSSHEET:
					{
						double	convfrom=0,convto=0,convinc=0;

						nconv = sscanf(pStartParam+1,"(%lf,%lf,%lf)",&convfrom,&convto,&convinc);
						if (nconv < LM_NPARAMS_MASSSHEET*3) {
							fprintf(stderr,"Parameter conversion failed for lens type %d. params were: <%s>\n",iCompType,pStartParam+1);
							iStatus = 1;
							break;
						}

						iStatus = lm_CreateLMComp_MassSheet(pLens,convfrom,convto,convinc);
					}
					break;

				case	LM_SIE:
					{
						double	xoff=0,yoff=0,critfrom=0,critto=0,critinc=0,axfrom=0,axto=0,axinc=0;
						double	angfrom=0,angto=0,anginc=0;

						nconv = sscanf(pStartParam+1,"(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf)",
                                       &xoff,&yoff,&critfrom,&critto, &critinc, &axfrom,&axto,&axinc,&angfrom,&angto,&anginc);

						if (nconv < LM_NPARAMS_SIE*3+2) {
							fprintf(stderr,"Parameter conversion failed for lens type %d. params were: <%s>\n",iCompType,pStartParam+1);
							iStatus = 1;
							break;
						}

						iStatus = lm_CreateLMComp_SIE(pLens,xoff,yoff,critfrom,critto,critinc,axfrom,axto,axinc, angfrom,angto,anginc);
					}
					break;

				case	LM_SPEMD:
					{
						double	xoff=0,yoff=0,critfrom=0,critto=0,critinc=0,ellfrom=0,ellto=0,ellinc=0;
						double	angfrom=0,angto=0,anginc=0,corefrom=0,coreto=0,coreinc=0,gammfrom=0,gammto=0,gamminc=0;
                        double  centerXfrom =0, centerXto=0, centerXinc = 0;
                        double  centerYfrom =0, centerYto=0, centerYinc = 0;

						nconv = sscanf(pStartParam+1,"(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf)",
                                       &xoff,&yoff,&critfrom,&critto, &critinc,&ellfrom,&ellto,&ellinc,
                                       &angfrom,&angto,&anginc,&gammfrom,&gammto,&gamminc,&corefrom,&coreto,&coreinc, &centerXfrom, &centerXto, &centerXinc, &centerYfrom, &centerYto, &centerYinc);

						if (nconv < LM_NPARAMS_SPEMD*3+2) {
							fprintf(stderr,"Parameter conversion failed for lens type %d. params were: <%s>\n",iCompType,pStartParam+1);
							iStatus = 1;
							break;
						}

						iStatus = lm_CreateLMComp_SPEMD(pLens,xoff,yoff,critfrom,critto,critinc,ellfrom,ellto,ellinc,
							angfrom,angto,anginc,gammfrom,gammto,gamminc,corefrom,coreto,coreinc, centerXfrom, centerXto, centerXinc, centerYfrom, centerYto, centerYinc);
					}
					break;

				case	LM_PTMASS:
					{
						double	xoff=0,yoff=0,critfrom=0,critto=0,critinc=0;

						nconv = sscanf(pStartParam+1,"(%lf,%lf,%lf,%lf,%lf)",&xoff,&yoff,&critfrom,&critto, &critinc);

						if (nconv < LM_NPARAMS_PTMASS*3 + 2) {
							fprintf(stderr,"Parameter conversion failed for lens type %d. params were: <%s>\n",iCompType,pStartParam+1);
							iStatus = 1;
							break;
						}
						iStatus = lm_CreateLMComp_PtMass(pLens,xoff,yoff,critfrom,critto,critinc);
					}
					break;

				case	LM_EXPDISC:
					{
						double xoff=0,yoff=0,kapfrom=0,kapto=0,kapinc=0,incfrom=0,incto=0,incinc=0;
						double angfrom=0,angto=0,anginc=0,scalefrom=0,scaleto=0,scaleinc=0;

						nconv = sscanf(pStartParam+1,"(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf)",
                                       &xoff,&yoff,&kapfrom,&kapto,&kapinc,&incfrom,&incto,&incinc,
                                       &angfrom,&angto,&anginc,&scalefrom,&scaleto,&scaleinc);
						if (nconv < LM_NPARAMS_EXPDISC*3 + 2) {
							fprintf(stderr,"Parameter conversion failed for lens type %d. params were: <%s>\n",iCompType,pStartParam+1);
							iStatus = 1;
							break;
						}
						iStatus = lm_CreateLMComp_ExpDisc(pLens,xoff,yoff,kapfrom,kapto,kapinc,incfrom,incto,incinc,
							angfrom,angto,anginc,scalefrom,scaleto,scaleinc);
					}
					break;
				case	LM_ISO_SPHERE:
					{
						double	xoff=0,yoff=0,critfrom=0,critto=0,critinc=0;

						nconv = sscanf(pStartParam+1,"(%lf,%lf,%lf,%lf,%lf)",&xoff,&yoff,&critfrom,&critto,&critinc);
						if (nconv < LM_NPARAMS_SIS*3 +2) {
							fprintf(stderr,"Parameter conversion failed for lens type %d. params were: <%s>\n",iCompType,pStartParam+1);
							iStatus = 1;
							break;
						}
						iStatus = lm_CreateLMComp_SIS(pLens,xoff,yoff,critfrom,critto,critinc);
					}
					break;
				case	LM_NFW:
					{
						double	xoff=0,yoff=0,massfrom=0,massto=0,massinc=0,scalefrom=0,scaleto=0,scaleinc=0;
						double	ellfrom=0,ellto=0,ellinc=0, angfrom=0,angto=0,anginc=0;

						nconv = sscanf(pStartParam+1,"(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf)",
                                       &xoff,&yoff,&massfrom,&massto,&massinc,&scalefrom,&scaleto,&scaleinc,
                                       &ellfrom,&ellto,&ellinc,&angfrom,&angto,&anginc);
						if (nconv < LM_NPARAMS_NFW*3+2) {
							fprintf(stderr,"Parameter conversion failed for lens type %d. params were: <%s>\n",iCompType,pStartParam+1);
							iStatus = 1;
							break;
						}
						iStatus = lm_CreateLMComp_NFW(pLens,xoff,yoff,massfrom,massto,massinc,scalefrom,scaleto,scaleinc,
							ellfrom,ellto,ellinc,angfrom,angto,anginc);
					}
					break;
				case	LM_SERSIC:
					{
						double xoff=0,yoff=0,kapfrom=0,kapto=0,kapinc=0,incfrom=0,incto=0,incinc=0;
						double angfrom=0,angto=0,anginc=0,scalefrom=0,scaleto=0,scaleinc=0,mfrom=0,mto=0,minc=0;

						nconv = sscanf(pStartParam+1,"(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf)",
                                       &xoff,&yoff,&kapfrom,&kapto,&kapinc,&incfrom,&incto,&incinc,
                                       &angfrom,&angto,&anginc,&scalefrom,&scaleto,&scaleinc,&mfrom,&mto,&minc);
						if (nconv < LM_NPARAMS_SERSIC*3 + 2) {
							fprintf(stderr,"Parameter conversion failed for lens type %d. params were: <%s>\n",iCompType,pStartParam+1);
							iStatus = 1;
							break;
						}
						iStatus = lm_CreateLMComp_Sersic(pLens,xoff,yoff,kapfrom,kapto,kapinc,incfrom,incto,incinc,
							angfrom,angto,anginc,scalefrom,scaleto,scaleinc,mfrom,mto,minc);
					}
					break;
				case	LM_FERRERS:
					{
						double	xoff=0,yoff=0,critfrom=0,critto=0,critinc=0,Afrom=0,Ato=0,Ainc=0;
						double	Bfrom=0,Bto=0,Binc=0,AngFrom=0,AngTo=0,AngInc=0,

						nconv = sscanf(pStartParam+1,"(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf)",
                                       &xoff,&yoff,&critfrom,&critto,&critinc,&Afrom,&Ato,&Ainc,&Bfrom,&Bto,&Binc,
                                       &AngFrom,&AngTo,&AngInc);
						if (nconv < LM_NPARAMS_FERRERS*3 +2) {
							fprintf(stderr,"Parameter conversion failed for lens type %d. params were: <%s>\n",iCompType,pStartParam+1);
							iStatus = 1;
							break;
						}
						iStatus = lm_CreateLMComp_Ferrers(pLens,xoff,yoff,critfrom,critto,critinc,Afrom,Ato,Ainc,Bfrom,Bto,Binc,AngFrom,AngTo,AngInc);
					}
					break;
				case	LM_USERDEF:
					{
						double	xoff=0,yoff=0,critfrom=0,critto=0,critinc=0;
						double	angfrom=0,angto=0,anginc=0,fPixAngSiz=0;
						char	strFilex[BUFF_SIZE]="",strFiley[BUFF_SIZE]="";
						lv_image_t	*pDeflx=NULL, *pDefly=NULL;

						nconv = sscanf(pStartParam+1,"(%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%s %s)",
                                       &xoff,&yoff,&critfrom,&critto,&critinc,&angfrom,&angto,&anginc,&fPixAngSiz,
                                       strFilex,strFiley);
						if (nconv < LM_NPARAMS_USERDEF*3 +5) {
							fprintf(stderr,"Parameter conversion failed for lens type %d. found %d params. params were: <%s>\n",iCompType,nconv,pStartParam+1);
							fprintf(stderr,"Filex: <%s>, filey: <%s>\n",strFilex,strFiley);
							iStatus = 1;
							break;
						}
						sprintf(strMessage,"Reading user-defined deflection fits file %s. Pixel size: %g",strFilex,fPixAngSiz);
						TRACE(LOG_HIGH_PRI,strMessage);
						pDeflx = lv_read_img_file(strFilex,fPixAngSiz);
						if (pDeflx == NULL) {
							sprintf("Unable to read file %s\n",strFilex);
							LOG_ERR(strMessage);
							iStatus=1;
							break;
						}

						sprintf(strMessage,"Reading user-defined deflection fits file %s. Pixel size: %g",strFiley,fPixAngSiz);
						TRACE(LOG_HIGH_PRI,strMessage);
						pDefly = lv_read_img_file(strFiley,fPixAngSiz);
						if (pDefly == NULL) {
							sprintf("Unable to read file %s\n",strFiley);
							LOG_ERR(strMessage);
							iStatus=1;
							break;
						}
						iStatus = lm_CreateLMComp_Userdef(pLens,xoff,yoff,critfrom,critto,critinc,angfrom,angto,anginc,pDeflx,pDefly);
					}
					break;
				default:
					fprintf(stderr,"lm_ReadLensFile: You forgot to handle lens type %d, you loser.\n",iCompType);
			}
			if (iStatus != 0) {
				LOG_ERR("Problem reading lens params. Aborting.");
				pLens = NULL;
				goto EXIT;
			}
		}
		else {
/*
			printf("Skipping line: %s",buffer);
*/
		}
	}

EXIT:
	if (fp !=NULL) {
		fclose(fp);
	}
	TRACE_OUT;
	return pLens;
}

/***************************
Function:
Description:
Arguments:
Returns:
****************************/
int	lm_decodeLensCompType(char *strName) {
	if (strcmp(strName,LM_NAME_PTMASS)==0) {
		return LM_PTMASS;
	}
	if (strcmp(strName,LM_NAME_PIEP)==0) {
		return LM_PIEP;
	}
	if (strcmp(strName,LM_NAME_SPEMD)==0) {
		return LM_SPEMD;
	}
	if (strcmp(strName,LM_NAME_NFW)==0) {
		return LM_NFW;
	}
	if (strcmp(strName,LM_NAME_SIE)==0) {
		return LM_SIE;
	}
	if (strcmp(strName,LM_NAME_EXTSH)==0) {
		return LM_EXTSHEAR;
	}
	if (strcmp(strName,LM_NAME_MASSSHT)==0) {
		return LM_MASSSHEET;
	}
	if (strcmp(strName,LM_NAME_SIS)==0) {
		return LM_ISO_SPHERE;
	}
	if (strcmp(strName,LM_NAME_EXPDISC)==0) {
		return LM_EXPDISC;
	}
	if (strcmp(strName,LM_NAME_SERSIC)==0) {
		return LM_SERSIC;
	}
	if (strcmp(strName,LM_NAME_FERRERS)==0) {
		return LM_FERRERS;
	}
	if (strcmp(strName,LM_NAME_USERDEF)==0) {
		return LM_USERDEF;
	}
	return 0;
}


/***************************
Function:	lm_CreateLMComp_Userdef
Description:	Create a lens model component with user defined deflection angles
Arguments:
Returns:
****************************/
static int	lm_CreateLMComp_Userdef(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t critfrom, real_t critto, real_t critinc, real_t AFrom, real_t ATo, real_t AInc, lv_image_t *pDefX, lv_image_t *pDefY) {
	int	iStatus =0;
	real_t	fromparams[LM_NPARAMS_FERRERS];
	real_t	toparams[LM_NPARAMS_FERRERS];
	real_t	incparams[LM_NPARAMS_FERRERS];
	char	*strNames[LM_NPARAMS_FERRERS] = {"Userdef_Kappa", "Userdef_Angle"};

	TRACE_IN(lm_CreateLMComp_Userdef);

	fromparams[0] = critfrom;
	toparams[0] = critto;
	incparams[0] = critinc;
	fromparams[1] = AFrom;
	toparams[1] = ATo;
	incparams[1] = AInc;

	pLens->iNumParameters +=LM_NPARAMS_USERDEF+2;
	iStatus = lm_InitLensModel(pLens, LM_USERDEF, fXoffset, fYoffset, LM_NPARAMS_USERDEF, fromparams, toparams, incparams,strNames);

	(pLens->pComponent[g_iCompNum-1]).pDeflX = pDefX;
	(pLens->pComponent[g_iCompNum-1]).pDeflY = pDefY;

	sprintf(strMessage,"Created lens model type %d with parameters: offset (%g,%g). Kappa from %g, to %g, inc %g. Angle from %g, to %g inc %g",LM_USERDEF,fXoffset,fYoffset, critfrom,critto,critinc, AFrom, ATo, AInc);
	TRACE(LOG_HIGH_PRI,strMessage);

	TRACE_OUT;
	return iStatus;
}


/***************************
Function:	lm_CreateLMComp_Ferrers
Description:	Create a Ferrers bar lens model component
Arguments:
Returns:
****************************/
static int	lm_CreateLMComp_Ferrers(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t critfrom, real_t critto, real_t critinc, real_t AFrom, real_t ATo, real_t AInc, real_t BFrom, real_t BTo, real_t BInc, real_t Angfrom, real_t Angto, real_t Anginc) {
	int	iStatus =0;
	real_t	fromparams[LM_NPARAMS_FERRERS];
	real_t	toparams[LM_NPARAMS_FERRERS];
	real_t	incparams[LM_NPARAMS_FERRERS];
	char	*strNames[LM_NPARAMS_FERRERS] = {"Ferrers_Kappa", "Ferrers_Axis_A", "Ferrers_Axis_B", "Ferrers_OrientAngle"};

	TRACE_IN(lm_CreateLMComp_Ferrers);

	fromparams[0] = critfrom;
	toparams[0] = critto;
	incparams[0] = critinc;
	fromparams[1] = AFrom;
	toparams[1] = ATo;
	incparams[1] = AInc;
	fromparams[2] = BFrom;
	toparams[2] = BTo;
	incparams[2] = BInc;
	fromparams[3] = Angfrom;
	toparams[3] = Angto;
	incparams[3] = Anginc;

	pLens->iNumParameters +=LM_NPARAMS_FERRERS+2;
	iStatus = lm_InitLensModel(pLens, LM_FERRERS, fXoffset, fYoffset, LM_NPARAMS_FERRERS, fromparams, toparams, incparams,strNames);

	sprintf(strMessage,"Created lens model type %d with parameters: offset (%g,%g). Kappa from %g, to %g, inc %g. Axis A from %g, to %g, inc %g. Axis B from %g, to %g, inc %g. Angle from %g, to %g inc %g",LM_FERRERS,fXoffset,fYoffset, critfrom,critto,critinc, AFrom, ATo, AInc, BFrom, BTo, BInc, Angfrom, Angto, Anginc);
	TRACE(LOG_HIGH_PRI,strMessage);

	TRACE_OUT;
	return iStatus;
}

/***************************
Function:	lm_CreateLMComp_SIS
Description:	Create a singluar isothermal sphere lens model component
Arguments:
Returns:
****************************/
int	lm_CreateLMComp_SIS(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t critfrom, real_t critto, real_t critinc) {
	int	iStatus =0;
	real_t	fromparams[LM_NPARAMS_SIS];
	real_t	toparams[LM_NPARAMS_SIS];
	real_t	incparams[LM_NPARAMS_SIS];
	char	*strNames[LM_NPARAMS_SIS] = {"Critrad"};

	TRACE_IN(lm_CreateLMComp_SIS);

	fromparams[0] = critfrom;
	toparams[0] = critto;
	incparams[0] = critinc;

	pLens->iNumParameters +=LM_NPARAMS_SIS+2;

	iStatus = lm_InitLensModel(pLens, LM_ISO_SPHERE, fXoffset, fYoffset, LM_NPARAMS_SIS, fromparams, toparams, incparams,strNames);

	sprintf(strMessage,"Created lens model type %d with parameters: offset (%g,%g). Critrad from %g, to %g, inc %g.",LM_ISO_SPHERE,fXoffset,fYoffset, critfrom,critto,critinc);
	TRACE(LOG_HIGH_PRI,strMessage);

	TRACE_OUT;
	return iStatus;
}

 
/***************************
Function:
Description:
Arguments:
Returns:
****************************/
int	lm_CreateLMComp_ExtShear(lv_lensmodel_t *pLens, real_t fShearFrom, real_t fShearTo, real_t fShearInc,
	real_t fAnlgeFrom, real_t fAngleTo, real_t fAngleInc) {

	int		iStatus=0;
	real_t	fromparams[LM_NPARAMS_EXTSH];
	real_t	toparams[LM_NPARAMS_EXTSH];
	real_t	incparams[LM_NPARAMS_EXTSH];
	char	*strNames[LM_NPARAMS_EXTSH] = {"Shear","Sh_Angle"};

	TRACE_IN(lm_CreateLMComp_ExtShear);

	fromparams[0] = fShearFrom;
	fromparams[1] = fAnlgeFrom;
	toparams[0] = fShearTo;
	toparams[1] = fAngleTo;
	incparams[0] = fShearInc;
	incparams[1] = fAngleInc;

	if (fShearFrom == 0.0 && fShearFrom < fShearTo && fAnlgeFrom < fAngleTo) {
		fprintf(stderr,"lm_CreateLMComp_ExtShear: WARNING: zero shear used for mulitple shear angles. This is redundant\n");
	}

	pLens->iNumParameters +=LM_NPARAMS_EXTSH;

	iStatus = lm_InitLensModel(pLens, LM_EXTSHEAR, 0, 0, LM_NPARAMS_EXTSH, fromparams, toparams, incparams,strNames);
	sprintf(strMessage,"Created lens model type %d with params: sh from: %g, sh to: %g, sh inc: %g, angle from: %g, ang to: %g, ang inc: %g"
				,LM_EXTSHEAR,fShearFrom,fShearTo,fShearInc,fAnlgeFrom,fAngleTo,fAngleInc);
	TRACE(LOG_HIGH_PRI,strMessage);
	TRACE_OUT;
	return iStatus;
}


/***************************
Function:
Description:
Arguments:
Returns:
****************************/
int	lm_CreateLMComp_NFW(lv_lensmodel_t *pLens, real_t fXoffset, real_t fYoffset, real_t fMassScaleFrom,
		real_t fMassScaleTo, real_t fMassScaleInc, real_t fScaleLenFrom, real_t fScaleLenTo, real_t fScaleLenInc,
		real_t fEllipFrom, real_t fEllipTo, real_t fEllipInc, real_t fAngleFrom, real_t fAngleTo, real_t fAngleInc) {
	int		iStatus=0;
	real_t	fromparams[LM_NPARAMS_NFW];
	real_t	toparams[LM_NPARAMS_NFW];
	real_t	incparams[LM_NPARAMS_NFW];
	char	*strNames[LM_NPARAMS_NFW] = {"Mass_scale","Scale_length","Ellipticity","Angle"};

	TRACE_IN(lm_CreateLMComp_NFW);

	if (fEllipFrom < 0 || fEllipTo >= 1.0) {
		sprintf(strMessage,"ERROR: Ellipticity must be between 0 and 1. Values were: %f to %f",fEllipFrom,fEllipTo);
		TRACE(LOG_ERROR_PRI,strMessage);
		iStatus = 1;
		goto EXIT;
	}

	if (fEllipFrom == 0 && fAngleFrom < fAngleTo) {
		sprintf(strMessage,"lm_CreateLMComp_NFW: WARNING: nonzero orientation range for zero inclination angle. This is redundant.\n");
		fprintf(stderr,strMessage);
		TRACE(LOG_HIGH_PRI,strMessage);
	}

	fromparams[0] = fMassScaleFrom;
	fromparams[1] = fScaleLenFrom;
	fromparams[2] = fEllipFrom;
	fromparams[3] = fAngleFrom;
	toparams[0] = fMassScaleTo;
	toparams[1] = fScaleLenTo;
	toparams[2] = fEllipTo;
	toparams[3] = fAngleTo;
	incparams[0] = fMassScaleInc;
	incparams[1] = fScaleLenInc;
	incparams[2] = fEllipInc;
	incparams[3] = fAngleInc;

	pLens->iNumParameters +=LM_NPARAMS_NFW+2;

	iStatus = lm_InitLensModel(pLens, LM_NFW,fXoffset,fYoffset, LM_NPARAMS_NFW, fromparams, toparams, incparams,strNames);

	sprintf(strMessage,"Created lens model type %d with offset (%g,%g). params: mass from: %g, mass to: %g, mass inc: %g, scale_len from: %g, to: %g, inc: %g, ellipticity from: %g, to: %g inc: %g, angle from: %g, to: %g, inc: %g",LM_NFW,fXoffset,fYoffset,
		fMassScaleFrom,fMassScaleTo,fMassScaleInc,fScaleLenFrom,fScaleLenTo,fScaleLenInc,fEllipFrom,fEllipTo,fEllipInc,
		fAngleFrom,fAngleTo,fAngleInc);
	TRACE(LOG_HIGH_PRI,strMessage);

EXIT:
	TRACE_OUT;
	return iStatus;
}


/***************************
Function:
Description:
Arguments:
Returns:
****************************/
int	lm_CreateLMComp_PtMass(lv_lensmodel_t *pLens, real_t fXoffset, real_t fYoffset, real_t fMassFrom,
	real_t fMassTo, real_t fMassInc) {
	int		iStatus=0;
	char	*strNames[LM_NPARAMS_PTMASS] = {"Critrad"};

	TRACE_IN(lm_CreateLMComp_PtMass);

	pLens->iNumParameters += LM_NPARAMS_PTMASS+2;

	iStatus = lm_InitLensModel(pLens, LM_PTMASS, fXoffset, fYoffset, LM_NPARAMS_PTMASS, &fMassFrom, &fMassTo, &fMassInc,strNames);

	sprintf(strMessage,"Created lens model type %d with offset (%g,%g). params: mass from: %g, mass to: %g, mass inc: %g"
				,LM_PTMASS,fXoffset,fYoffset,fMassFrom,fMassTo,fMassInc);
	TRACE(LOG_HIGH_PRI,strMessage);

	TRACE_OUT;
	return iStatus;
}


/***************************
Function:
Description:
Arguments:
Returns:
****************************/
int	lm_CreateLMComp_MassSheet(lv_lensmodel_t *pLens, real_t fScaleFrom, real_t fScaleTo, real_t fScaleInc) {

	int		iStatus=0;
	char	*strNames[LM_NPARAMS_MASSSHEET] = {"Convergence"};

	pLens->iNumParameters +=1;

	TRACE_IN(lm_CreateLMComp_MassSheet);

	iStatus = lm_InitLensModel(pLens, LM_MASSSHEET, 0, 0, LM_NPARAMS_MASSSHEET, &fScaleFrom, &fScaleTo, &fScaleInc,strNames);

	TRACE_OUT;
	return iStatus;
}


/***************************
Function:
Description:
Arguments:
Returns:
****************************/
static int	lm_CreateLMComp_SPEMD(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t fMassScaleFrom, 
		real_t fMassScaleTo, real_t fMassScaleInc, real_t fEllipFrom, real_t fEllipTo, real_t fEllipInc,
		real_t fAngleFrom, real_t fAngleTo, real_t fAngleInc, real_t fGammaFrom, real_t fGammaTo,
		real_t fGammaInc, real_t fCoreFrom, real_t fCoreTo, real_t fCoreInc,
		real_t fCenterXFrom, real_t fCenterXTo, real_t fCenterXInc,
		real_t fCenterYFrom, real_t fCenterYTo, real_t fCenterYInc) {

	int		iStatus=0;
	real_t	fromparams[LM_NPARAMS_SPEMD];
	real_t	toparams[LM_NPARAMS_SPEMD];
	real_t	incparams[LM_NPARAMS_SPEMD];
	char	*strNames[LM_NPARAMS_SPEMD] = {"Critrad","Ellipticity","Orient_Angle","Gamma","Core_rad", "CenterX", "CenterY"};

	TRACE_IN(lm_CreateLMComp_SPEMD);

	if (fGammaFrom < 0.0 || fGammaTo > 1.0) {
		sprintf(strMessage,"ERROR: Gamma must be between 0 and 1. Values were: %f to %f",fGammaFrom,fGammaTo);
		TRACE(LOG_ERROR_PRI,strMessage);
		iStatus =1;
		goto EXIT;
	}

	fromparams[0] = fMassScaleFrom;
	fromparams[1] = fEllipFrom;
	fromparams[2] = fAngleFrom;
	fromparams[3] = fGammaFrom;
	fromparams[4] = fCoreFrom;
	fromparams[5] = fCenterXFrom; 
	fromparams[6] = fCenterYFrom; 

	toparams[0] = fMassScaleTo;
	toparams[1] = fEllipTo;
	toparams[2] = fAngleTo;
	toparams[3] = fGammaTo;
	toparams[4] = fCoreTo;
	toparams[5] = fCenterXTo; 
	toparams[6] = fCenterYTo; 

	incparams[0] = fMassScaleInc;
	incparams[1] = fEllipInc;
	incparams[2] = fAngleInc;
	incparams[3] = fGammaInc;
	incparams[4] = fCoreInc;
	incparams[5] = fCenterXInc; 
	incparams[6] = fCenterYInc; 

	pLens->iNumParameters +=LM_NPARAMS_SPEMD+2;

	iStatus = lm_InitLensModel(pLens, LM_SPEMD, fXoffset, fYoffset, LM_NPARAMS_SPEMD , fromparams, toparams, incparams,strNames); 

	sprintf(strMessage,"Created SPEMD lens model component. offset: (%.2f,%.2f), Critrad: %.3f-%.3f +%.3f, Ellip: %.3f-%.3f +%.3f, Angle: %.1f-%.1f +%.1f, Gamma: %.3f-%.3f +%.3f, Core: %.3f-%.3f +%.3f, CenterX: %.3f-%.3f +%.3f, CenterY: %.3f-%.3f +%.3f",fXoffset,fYoffset,fMassScaleFrom, fMassScaleTo, fMassScaleInc, fEllipFrom, fEllipTo, fEllipInc, fAngleFrom, fAngleTo, fAngleInc, fGammaFrom,fGammaTo,fGammaInc,fCoreFrom, fCoreTo, fCoreInc, fCenterXFrom, fCenterXTo, fCenterXInc,fCenterYFrom, fCenterYTo, fCenterYInc );
	TRACE(LOG_HIGH_PRI,strMessage);

	/*
	LOG_ERR("this lens model is not available");
	iStatus=1;
	*/


EXIT:
	TRACE_OUT;
	return iStatus;
}


/***************************
Function:
Description:
Arguments:	fMassScaleFrom: Kappa, fIncFrom: axis ratio of light dist,
			fAngleFrom: orientation angle for non-circ profile,
			fScaleLenFrom: scale len in arcseconds
Returns:
****************************/
int	lm_CreateLMComp_Sersic(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t fMassScaleFrom,
		real_t fMassScaleTo, real_t fMassScaleInc, real_t fIncFrom, real_t fIncTo, real_t fIncInc,
		real_t fAngleFrom, real_t fAngleTo, real_t fAngleInc, real_t fScaleLenFrom, real_t fScaleLenTo, real_t fScaleLenInc,
		real_t fMFrom, real_t fMTo, real_t fMInc) {

	int		iStatus=1;
	real_t	fromparams[LM_NPARAMS_SERSIC];
	real_t	toparams[LM_NPARAMS_SERSIC];
	real_t	incparams[LM_NPARAMS_SERSIC];
	char	*strNames[LM_NPARAMS_SERSIC] = {"Sersic_kappa","Sersic_Axratio","Sersic_OrientAngle","Sersic_scale","R_Power"};

	TRACE_IN(lm_CreateLMComp_Sersic);

	if (fIncFrom == 1 && fAngleFrom < fAngleTo) {
		fprintf(stderr,"lm_CreateLMComp_Sersic: WARNING: nonzero orientation range for zero ellipticity (axratio=1). This is redundant\n");
	}
	if (fIncFrom <= 0.0 || fIncTo > 1.0) {
		sprintf(strMessage,"ERROR: axis ratio parameter must be 0 < q <= 1.0. Given: %g to %g",fIncFrom,fIncTo);
		TRACE(LOG_ERROR_PRI,strMessage);
		goto EXIT;
	}
	if (fScaleLenFrom <= 0.0) {
		sprintf(strMessage,"ERROR: disc scale must be > 0. Given: %g",fScaleLenFrom);
		TRACE(LOG_ERROR_PRI,strMessage);
		goto EXIT;
	}

	fromparams[0] = fMassScaleFrom;
	fromparams[1] = fIncFrom;
	fromparams[2] = fAngleFrom;
	fromparams[3] = fScaleLenFrom;
	fromparams[4] = fMFrom;

	toparams[0] = fMassScaleTo;
	toparams[1] = fIncTo;
	toparams[2] = fAngleTo;
	toparams[3] = fScaleLenTo;
	toparams[4] = fMTo;

	incparams[0] = fMassScaleInc;
	incparams[1] = fIncInc;
	incparams[2] = fAngleInc;
	incparams[3] = fScaleLenInc;
	incparams[4] = fMInc;

	pLens->iNumParameters +=LM_NPARAMS_SERSIC+2;

	iStatus = lm_InitLensModel(pLens, LM_SERSIC, fXoffset, fYoffset, LM_NPARAMS_SERSIC, fromparams, toparams, incparams,strNames);

	if (iStatus == 0) {
		sprintf(strMessage,"Created lens model type %d with params: offset(%g,%g), kap from: %g, to: %g, inc: %g, axratio from: %g, to: %g, inc: %g. Angle from: %g, to %g, inc %g. Scale len from: %g, to %g, inc %g. M from: %g, to: %g, inc: %g" ,LM_SERSIC,fXoffset,fYoffset,fMassScaleFrom,fMassScaleTo,fMassScaleInc,fIncFrom,fIncTo,fIncInc, fAngleFrom,fAngleTo,fAngleInc,fScaleLenFrom,fScaleLenTo,fScaleLenInc,fMFrom,fMTo,fMInc);
		TRACE(LOG_HIGH_PRI,strMessage);
	}

EXIT:
	TRACE_OUT;
	return iStatus;
}


/***************************
Function:
Description:
Arguments:
Returns:
****************************/
int	lm_CreateLMComp_ExpDisc(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t fMassScaleFrom,
		real_t fMassScaleTo, real_t fMassScaleInc, real_t fIncFrom, real_t fIncTo, real_t fIncInc,
		real_t fAngleFrom, real_t fAngleTo, real_t fAngleInc, real_t fScaleLenFrom, real_t fScaleLenTo, real_t fScaleLenInc) {

	int		iStatus=1;
	real_t	fromparams[LM_NPARAMS_EXPDISC];
	real_t	toparams[LM_NPARAMS_EXPDISC];
	real_t	incparams[LM_NPARAMS_EXPDISC];
	char	*strNames[LM_NPARAMS_EXPDISC] = {"Disc_kappa","Incl_axratio","Orient_Angle","Disc_scale"};

	TRACE_IN(lm_CreateLMComp_ExpDisc);

	if (fIncFrom == 0 && fAngleFrom < fAngleTo) {
		fprintf(stderr,"lm_CreateLMComp_ExpDisc: WARNING: nonzero orientation range for zero inclination angle. This is redundant\n");
	}
	if (fIncFrom < 0.0 || fIncTo >= 1.0) {
		sprintf(strMessage,"ERROR: inclination parameter must be 0 <= q < 1.0. Given: %g to %g",fIncFrom,fIncTo);
		TRACE(LOG_ERROR_PRI,strMessage);
		goto EXIT;
	}
	if (fScaleLenFrom <= 0.0) {
		sprintf(strMessage,"ERROR: disc scale must be > 0. Given: %g",fScaleLenFrom);
		TRACE(LOG_ERROR_PRI,strMessage);
		goto EXIT;
	}

	fromparams[0] = fMassScaleFrom;
	fromparams[1] = fIncFrom;
	fromparams[2] = fAngleFrom;
	fromparams[3] = fScaleLenFrom;

	toparams[0] = fMassScaleTo;
	toparams[1] = fIncTo;
	toparams[2] = fAngleTo;
	toparams[3] = fScaleLenTo;

	incparams[0] = fMassScaleInc;
	incparams[1] = fIncInc;
	incparams[2] = fAngleInc;
	incparams[3] = fScaleLenInc;

	pLens->iNumParameters +=LM_NPARAMS_EXPDISC+2;

	iStatus = lm_InitLensModel(pLens, LM_EXPDISC, fXoffset, fYoffset, LM_NPARAMS_EXPDISC, fromparams, toparams, incparams,strNames);

	if (iStatus == 0) {
		sprintf(strMessage,"Created lens model type %d with params: kap from: %g, to: %g, inc: %g, axratio from: %g, to: %g, inc: %g. Angle from: %g, to %g, inc %g. Scale len from: %g, to %g, inc %g" ,LM_EXPDISC,fMassScaleFrom,fMassScaleTo,fMassScaleInc,fIncFrom,fIncTo,fIncInc, fAngleFrom,fAngleTo,fAngleInc,fScaleLenFrom,fScaleLenTo,fScaleLenInc);
		TRACE(LOG_HIGH_PRI,strMessage);
	}

EXIT:
	TRACE_OUT;
	return iStatus;
}

/***************************
Function:
Description:
Arguments:
Returns:
****************************/
int	lm_CreateLMComp_SIE(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t fMassScaleFrom, 
		real_t fMassScaleTo, real_t fMassScaleInc, real_t fAxFrom, real_t fAxTo, real_t fAxInc,
		real_t fAngleFrom, real_t fAngleTo, real_t fAngleInc) {

	int		iStatus=0;
	real_t	fromparams[LM_NPARAMS_SIE];
	real_t	toparams[LM_NPARAMS_SIE];
	real_t	incparams[LM_NPARAMS_SIE];
	char	*strNames[LM_NPARAMS_SIE] = {"Critrad","Axis_ratio","Orient_Angle"};

	TRACE_IN(lm_CreateLMComp_SIE);

	if (fAxFrom <= 0 || fAxTo > 1.0) {
		sprintf(strMessage,"ERROR: Axis_ratio must be > 0 and <= 1. Values were: %f to %f",fAxFrom,fAxTo);
		TRACE(LOG_ERROR_PRI,strMessage);
		iStatus = 1;
		goto EXIT;
	}

	fromparams[0] = fMassScaleFrom;
	fromparams[1] = fAxFrom;
	fromparams[2] = fAngleFrom;

	toparams[0] = fMassScaleTo;
	toparams[1] = fAxTo;
	toparams[2] = fAngleTo;

	incparams[0] = fMassScaleInc;
	incparams[1] = fAxInc;
	incparams[2] = fAngleInc;

	pLens->iNumParameters +=LM_NPARAMS_SIE+2;

	iStatus = lm_InitLensModel(pLens, LM_SIE, fXoffset, fYoffset, LM_NPARAMS_SIE, fromparams, toparams, incparams,strNames);
	if (iStatus==0) {
		sprintf(strMessage,"Created SIE lens model component. offset: (%.2f,%.2f), Critrad: %g-%g +%g Axratio: %g-%g +%g, Angle: %.1f-%.1f +%.1f",fXoffset,fYoffset,fMassScaleFrom, fMassScaleTo, fMassScaleInc, fAxFrom, fAxTo, fAxInc, fAngleFrom, fAngleTo, fAngleInc);
		TRACE(LOG_HIGH_PRI,strMessage);
	}

EXIT:
	TRACE_OUT;
	return iStatus;
}


/***************************
Function:
Description:
Arguments:
Returns:
****************************/
int	lm_CreateLMComp_PIEP(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t fMassScaleFrom,
	real_t fMassScaleTo, real_t fMassScaleInc, real_t fEllipFrom, real_t fEllipTo, real_t fEllipInc,
	real_t fAngleFrom, real_t fAngleTo, real_t fAngleInc, real_t fCoreFrom, real_t fCoreTo, real_t fCoreInc) {

	int		iStatus=0;
	real_t	fromparams[LM_NPARAMS_PIEP];
	real_t	toparams[LM_NPARAMS_PIEP];
	real_t	incparams[LM_NPARAMS_PIEP];
	char	*strNames[LM_NPARAMS_PIEP] = {"Critrad","Ellipticity","Orient_Angle","Core_rad"};

	TRACE_IN(lm_CreateLMComp_PIEP);
	if (fEllipFrom < 0 || fEllipTo >= 1.0) {
		sprintf(strMessage,"ERROR: Ellipticity must be between 0 and 1. Values were: %f to %f",fEllipFrom,fEllipTo);
		TRACE(LOG_ERROR_PRI,strMessage);
		iStatus = 1;
		goto EXIT;
	}

	if (fEllipFrom == 0 && fAngleFrom < fAngleTo) {
		fprintf(stderr,"lm_CreateLMComp_PIEP: WARNING: ellipticity is zero for a range of angles. This is redundant.\n");
	}

	fromparams[0] = fMassScaleFrom;
	fromparams[1] = fEllipFrom;
	fromparams[2] = fAngleFrom;
	fromparams[3] = fCoreFrom;
	toparams[0] = fMassScaleTo;
	toparams[1] = fEllipTo;
	toparams[2] = fAngleTo;
	toparams[3] = fCoreTo;
	incparams[0] = fMassScaleInc;
	incparams[1] = fEllipInc;
	incparams[2] = fAngleInc;
	incparams[3] = fCoreInc;

	pLens->iNumParameters +=LM_NPARAMS_PIEP+2;

	iStatus = lm_InitLensModel(pLens, LM_PIEP, fXoffset, fYoffset, LM_NPARAMS_PIEP, fromparams, toparams, incparams,strNames);

	sprintf(strMessage,"Created PIEP lens model component. offset: (%.2f,%.2f), Critrad: %g-%g +%g, Ellip: %.3f-%.3f +%.3f, Angle: %.1f-%.1f +%.1f, Core: %.3f-%.3f +%.3f",fXoffset,fYoffset,fMassScaleFrom, fMassScaleTo, fMassScaleInc, fEllipFrom, fEllipTo, fEllipInc, fAngleFrom, fAngleTo, fAngleInc, fCoreFrom, fCoreTo, fCoreInc);
	TRACE(LOG_HIGH_PRI,strMessage);

EXIT:
	TRACE_OUT;
	return iStatus;
}


/***************************
Function:	lm_InitLensModel
Description:	Initialise a lens model component with relevant parameters
Arguments:	pLens:	pointer to lens which was generated by lm_CreateLensModel or equivalent
			fXoffset, fYoffset: offset of centre of lens from centre of image plane
			rest: relevant lens parameters. fShearAngle in radians. fGalRadius not used.
			fCritRadius	in arcsec
Returns:	0:	success
			other: failure
****************************/
int	lm_InitLensModel(lv_lensmodel_t *pLens, int	iType, real_t fXoffset, real_t fYoffset, int iNumparams, real_t fParamsfrom[], real_t fParamsto[], real_t fParamsinc[],char *strNames[]) {

	int		iStatus=0,i;

	TRACE_IN(lv_InitLensModel);

	if (g_iCompNum == pLens->iNumComponents) {
		sprintf(strMessage,"Tried to create too many lens compoents. Current: %d. Max: %d",g_iCompNum+1,pLens->iNumComponents);
		TRACE(LOG_ERROR_PRI,strMessage);
		iStatus=1;
		goto EXIT;
	}

	if (iNumparams > MAX_LENSCOMP_PARAMS) {
		sprintf(strMessage,"Too many parameters in lens component. Given: %d, Max: %d",iNumparams,MAX_LENSCOMP_PARAMS);
		TRACE(LOG_ERROR_PRI,strMessage);
		iStatus=1;
		goto EXIT;
	}

	pLens->pComponent[g_iCompNum].iType = iType;
	pLens->pComponent[g_iCompNum].iNumParams = iNumparams;
	pLens->pComponent[g_iCompNum].fXoffset = fXoffset;
	pLens->pComponent[g_iCompNum].fYoffset = fYoffset;

	for (i=0; i< iNumparams; i++) {
		pLens->pComponent[g_iCompNum].fParameter[i] = fParamsfrom[i];
		pLens->pComponent[g_iCompNum].fParamfrom[i] = fParamsfrom[i];
		pLens->pComponent[g_iCompNum].fParamto[i] = fParamsto[i];
		pLens->pComponent[g_iCompNum].fParaminc[i] = fParamsinc[i];
		pLens->pComponent[g_iCompNum].strParamName[i] = strdup(strNames[i]);

		if (fParamsfrom[i] > fParamsto[i] ) {
			sprintf(strMessage,"ERROR: lens model type %d has from param (%g) greater than to param (%g) for component %d. Typo perhaps?\n",iType,fParamsfrom[i],fParamsto[i],g_iCompNum+1);
			fprintf(stderr,strMessage);
			LOG_ERR(strMessage)
			iStatus=1;
			goto EXIT;
		}
	}

	g_iCompNum++;

	sprintf(strMessage,"Created Lens Model type %d. This is %d out of %d",iType,g_iCompNum,pLens->iNumComponents);
	TRACE(LOG_MED_PRI, strMessage);

EXIT:
	TRACE_OUT;
	return iStatus;
}

/***************************
Function:	lm_CreateLensModel
Description:	allocate and initialise memory for a lens model structure
Arguments:	iNumComponents:	number of componets which make up total lens
Returns:	NULL: failure
			non-null: pointer to lv_lensmodel_t structure
****************************/
lv_lensmodel_t	*lm_CreateLensModel(int	iNumComp) {
	lv_lensmodel_t	*pLensModel = NULL;

	TRACE_IN(lm_CreateLensModel);
	pLensModel = (lv_lensmodel_t *) malloc(sizeof(lv_lensmodel_t));
	if (pLensModel == NULL) {
		fprintf(stderr,"lm_CreateLensModel: no malloc\n");
		goto EXIT;
	}

	/* init stuff */
	pLensModel->iNumComponents = iNumComp;
	pLensModel->iNumParameters =0;

	pLensModel->pComponent = calloc(iNumComp,sizeof(lv_lenscomp));
	if (pLensModel->pComponent == NULL) {
		fprintf(stderr,"lm_CreateLensModel: no malloc\n");
		free(pLensModel);
		pLensModel = NULL;
	}

EXIT:
	TRACE_OUT;
	return pLensModel;
}

/***************************
Function:	lm_FreeLensModel
Description:	free memory associated with the lens model
Arguments:	pLensModel:	 pointer to lens
Returns:	void
****************************/
void	lm_FreeLensModel(lv_lensmodel_t *pLensModel) {

	int i,j;
	lv_lenscomp *pComp;

	TRACE_IN(lv_FreeLensModel);
	if (pLensModel == NULL) {
		goto EXIT;
	}
	for (i=0; i< pLensModel->iNumComponents; i++) {
		pComp = pLensModel->pComponent+i;
		if (pComp->pDeflX != NULL) lv_free_image_struct(pComp->pDeflX);
		if (pComp->pDeflY != NULL) lv_free_image_struct(pComp->pDeflY);
		for (j=0; j< pComp->iNumParams; j++) {
			if (pComp->strParamName[j] != NULL) free(pComp->strParamName[j]);
		}
	}
	free(pLensModel->pComponent);
	free(pLensModel);

EXIT:
	TRACE_OUT;
}

/***************************
Function:	lm_nfw_mass
Description:	Calculate the analytical mass inside an NFW halo at radius x
Arguments:	x: radial distance in units of scale radius
Returns:	mass or 0 on error.
****************************/
static	double	lm_nfw_mass(double x) {
	double	logx_2=0,result=0;

	if (x < 0) {
		fprintf(stderr,"lm_nfw_mass: invalid x: %g. Must be > 0\n",x);
		return result;
	}

	if (x==0) return 0;
	logx_2 = log(0.5 * x);
	if (x < 1.0) {
		result = lm_arccosh(1/x)/sqrt(1 - x*x);
/*
		result = 2/sqrt(1 - x*x) * lm_arctanh(sqrt((1-x)/(1+x)));
*/
	}
	else if (x == 1.0) {
		result =  1.0;
	}
	else {	 /* x > 1 */
		result =  acos(1.0/x)/sqrt(x*x -1.);
	}

	return 4.0*(logx_2 + result);
}

/***************************
Function:	lm_arctanh
Description:	Compute the arc hyperbolic tangent analytically
Arguments:	x: input between 0 and 1.
Returns:	arcthan(x) or 0 on error
****************************/
static	double	lm_arctanh(double x) {
	if (x < -1 || x > 1.) {
		fprintf(stderr,"lm_arctanh: invalid x: %g. Must be 0 <= x <= 1\n",x);
		return 0;
	}
	return log(sqrt((1.+x)/(1.-x)));
}

/***************************
Function:	lm_arccosh
Description:	Compute the arc hyperbolic cosine analytically
Arguments:	x: input >= 1.
Returns:	arccosh(x) or 0 on error
****************************/
static	double	lm_arccosh(double x) {
	if (x < 1.) {
		fprintf(stderr,"lm_arccosh: invalid x: %g. Must be >= 1.\n",x);
		return 0;
	}
	return log(sqrt(x*x-1.)+x);
}

/***************************
Function:	lm_expdiscdeflx
Description:	calculate the x deflection for an exponential disc
Arguments:	x,y: coords in image plane. q: axis ratio of disc equal to cos(i)
			where i is the inclination angle. i=0 is face on
			global variables for scale len and axisratio are necessary for the
			numerical integration
Returns:	deflection angle (arcsec)
****************************/
real_t	lm_expdiscdeflx(real_t x, real_t y, real_t fAxRatio, real_t fScaleLen){

	static	gsl_integration_workspace   *workspace=NULL;
	gsl_function    gfunc;
	int	iResult=0;
	double	res=0,abserr=0;

	TRACE_IN(lm_expdiscdeflx);

	if (x==0) {
		goto EXIT;
	}
	gfunc.function = &lm_expdiscintx_gsl;
	gfunc.params = NULL;
	if (workspace==NULL) {
		workspace = gsl_integration_workspace_alloc(1000);
	}

	g_axratio = fAxRatio;
	g_scale = fScaleLen;
	g_tempx = x;
	g_tempy = y;

	iResult = gsl_integration_qag(&gfunc,0.0,1.0,0.0,1e-6,1000,GSL_INTEG_GAUSS61,workspace,&res,&abserr);
	if (iResult != 0) {
		sprintf(strMessage,"ERROR: integration problems. return val: %d, x: %g, y:%g, axratio: %g, scale: %g, res: %g, abserr: %g\n",
			iResult,x,y,fAxRatio,fScaleLen,res,abserr);
		LOG_ERR(strMessage);
	}
EXIT:
	TRACE_OUT;
	return (real_t) res*x;

/*
printf("calc x defl for (x,y) = (%g,%g). axratio: %g, scale %g\n",x,y,fAxRatio,fScaleLen);
printf("discintx0.5 is %g\n",lm_expdiscintx(0.5));
*/
/*
	return x*(qromo(lm_expdiscintx,0.0,0.1,midpnt) + qromo(lm_expdiscintx,0.1,1.0,midpnt));
*/
}

/***************************
Function:	lm_expdiscdefly
Description:	see lm_expdiscdeflx
Arguments:
Returns:
****************************/
real_t	lm_expdiscdefly(real_t x, real_t y, real_t fAxRatio, real_t fScaleLen){

	static	gsl_integration_workspace   *workspace=NULL;
	gsl_function    gfunc;
	int	iResult=0;
	double	res=0,abserr=0;

	TRACE_IN(lm_expdiscdefly)

	if (y==0) {
		return 0;
	}

	gfunc.function = &lm_expdiscinty_gsl;
	gfunc.params = NULL;
	if (workspace==NULL) {
		workspace = gsl_integration_workspace_alloc(1000);
	}

	g_axratio =fAxRatio;
	g_scale = fScaleLen;
	g_tempx = x;
	g_tempy = y;

	iResult = gsl_integration_qag(&gfunc,0.0,1.0,0.0,1e-6,1000,GSL_INTEG_GAUSS61,workspace,&res,&abserr);
	if (iResult != 0) {
		sprintf(strMessage,"ERROR: integration problems. return val: %d, x: %g, y:%g, axratio: %g, scale: %g, res: %g, abserr: %g\n",
			iResult,x,y,fAxRatio,fScaleLen,res,abserr);
		LOG_ERR(strMessage);
	}

	TRACE_OUT;
	return (real_t) res*y;
/*
	return y*(qromo(lm_expdiscinty,0.0,0.1,midpnt) + qromo(lm_expdiscinty,0.1,1.0,midpnt));
*/
}

/***************************
Function:	lm_expdiscintx
Description:	Integration function for qsimp for the deflection of an
				exponential disk
Arguments:	u: funtion variable
Returns:	function value
****************************/
/*
static	real_t	lm_expdiscintx(real_t u) {
	real_t	num=0,denom=0;

	denom = 1.0 - (1.0 - g_axratio*g_axratio)*u;
	num = u*(g_tempx*g_tempx + (g_tempy*g_tempy)/denom);
	return exp(-sqrt(num)/g_scale)/sqrt(denom);
}
*/
static	double	lm_expdiscintx_gsl(double u, void *params) {
	real_t	num=0,denom=0;

	denom = 1.0 - (1.0 - g_axratio*g_axratio)*u;
	num = u*(g_tempx*g_tempx + (g_tempy*g_tempy)/denom);
	return exp(-sqrt(num)/g_scale)/sqrt(denom);
}

/***************************
Function:	lm_expdiscinty
Description:	Integration function for qsimp for the deflection of an exponential disk
Arguments:	u: funtion variable
Returns:	function value
****************************/
/*
static	real_t	lm_expdiscinty(real_t u) {
	real_t	num=0,denom=0;

	denom = 1.0 - (1.0 - g_axratio*g_axratio)*u;
	num = u*(g_tempx*g_tempx + (g_tempy*g_tempy)/denom);
	return exp(-sqrt(num)/g_scale)/pow(denom,1.5);
}
*/
static	double	lm_expdiscinty_gsl(double u, void *params) {
	real_t	num=0,denom=0;

	denom = 1.0 - (1.0 - g_axratio*g_axratio)*u;
	num = u*(g_tempx*g_tempx + (g_tempy*g_tempy)/denom);
	return exp(-sqrt(num)/g_scale)/pow(denom,1.5);
}


/***************************
Function:
Description:
Arguments:
Returns:
****************************/
static	int	lm_deflCacheLookup(short iType, lmDeflCache *pCacheList, real_t x, real_t y, real_t fAxratio, real_t fScale, real_t fM, real_t *pDeflx, real_t *pDefly) {

	int	iRes = 0,iFound = FALSE,i,j;
	int	iDeltaX=0,iDeltaY=0;
	lmDeflCache	*pTemp=NULL,*pCache=NULL;
	real_t	fTempx, fTempy,t,u,xgrid,ygrid,fResx,fResy;
	int		x1,x2,y1,y2;

	TRACE_IN(lm_deflCacheLookup);

	/* find the cache for the current parameters */
	pTemp = pCacheList;
	if ( pTemp->pNext != NULL) {
		do {
			pTemp = pTemp->pNext;
			if(pTemp->param[0]==fAxratio && pTemp->param[1]==fScale && pTemp->param[2]==fM) {
				iFound=TRUE;
				pCache = pTemp;
				break;
			}
		} while (pTemp->pNext != NULL);
	}

	/* if there wasn't a cache for the current params, then we must make one */
	if (iFound==FALSE) {
		sprintf(strMessage,"Creating new cache for mass dist %d with params axratio %g and scale %g",iType,fAxratio,fScale);
		TRACE(LOG_HIGH_PRI,strMessage);

		/* first call to new cache params. Set things up */
		pCache = (lmDeflCache *) malloc(sizeof(lmDeflCache));
		if (pCache == NULL) {
			LOG_ERR("No malloc for new cache struct");
			iRes = 1;
			goto EXIT;
		}
		pTemp->pNext = pCache;
		pCache->pNext = NULL;
		pCache->iType = iType;
		pCache->pix_size = g_PixelResn;
		pCache->param[0] = fAxratio;
		pCache->param[1] = fScale;
		pCache->param[2] = fM;
		/* make the cache the size of a grid rotated 45 deg plus some extra space */
		pCache->dimension[0] = g_iSzImgx*1.41+20;
		pCache->dimension[1] = g_iSzImgy*1.41+20;
		/* make the cache an odd sized array so that zero is in the middle */
		if (pCache->dimension[0]%2==0) pCache->dimension[0] +=1;
		if (pCache->dimension[1]%2==0) pCache->dimension[1] +=1;
		/* make space for the x,y ordinate values */
		pCache->pValsX = calloc(pCache->dimension[0],sizeof(real_t));
		pCache->pValsY = calloc(pCache->dimension[1],sizeof(real_t));
		for (i=0; i<pCache->dimension[0]; i++)	pCache->pValsX[i] = (i-pCache->dimension[0]/2)*g_PixelResn;
		for (i=0; i<pCache->dimension[1]; i++)	pCache->pValsY[i] = (i-pCache->dimension[1]/2)*g_PixelResn;
		/* make space for the x,y deflection values */
		pCache->pDeflX = calloc(pCache->dimension[0]*pCache->dimension[1],sizeof(real_t));
		pCache->pDeflY = calloc(pCache->dimension[0]*pCache->dimension[1],sizeof(real_t));
		if (pCache->pValsX == NULL || pCache->pValsY ==NULL || pCache->pDeflX==NULL || pCache->pDeflY==NULL){
			LOG_ERR("No malloc");
			iRes = -1;
			goto EXIT;
		}
		/* fill up the cache! */
		sprintf(strMessage,"Filling %dx%d cache.",(int)pCache->dimension[0],(int)pCache->dimension[1]);
		TRACE(LOG_MED_PRI,strMessage);
		for (j=0; j<=pCache->dimension[1]/2; j++) {

			fTempy=(j-pCache->dimension[1]/2)*pCache->pix_size;

			for (i=0; i<=pCache->dimension[0]/2; i++) {

				fTempx = (i-pCache->dimension[0]/2)*pCache->pix_size;

				switch (iType) {
					case	LM_SERSIC:
						iRes = lm_CalcSersicDefl(fTempx,fTempy,fAxratio,fScale,fM,&fResx,&fResy);
						break;
					case	LM_EXPDISC:
						fResx = lm_expdiscdeflx(fTempx,fTempy,fAxratio,fScale);
						fResy = lm_expdiscdefly(fTempx,fTempy,fAxratio,fScale);
						break;
					default:
						sprintf(strMessage,"ERROR: Unknown lens type %d",iType);
						LOG_ERR(strMessage)
						iRes = 1;
						goto EXIT;
						break;
				}
				/* assume that the deflections are symmetric around both axes */
				pCache->pDeflX[j*(pCache->dimension[0]) + i] = fResx;
				pCache->pDeflX[j*(pCache->dimension[0]) + (pCache->dimension[0]-i-1)] = -fResx;
				pCache->pDeflX[(pCache->dimension[1]-j-1)*(pCache->dimension[0]) + (pCache->dimension[0]-i-1)] = -fResx;
				pCache->pDeflX[(pCache->dimension[1]-j-1)*(pCache->dimension[0]) + i] = fResx;
				pCache->pDeflY[j*(pCache->dimension[0]) + i] = fResy;
				pCache->pDeflY[j*(pCache->dimension[0]) + (pCache->dimension[0]-i-1)] = fResy;
				pCache->pDeflY[(pCache->dimension[1]-j-1)*(pCache->dimension[0]) + (pCache->dimension[0]-i-1)] = -fResy;
				pCache->pDeflY[(pCache->dimension[1]-j-1)*(pCache->dimension[0]) + i] = -fResy;
			}
		}
		TRACE(LOG_MED_PRI,"Done filling cache.");
/*
		img_DumpImage(pCache->pDeflX,pCache->dimension,"deflX");
		img_DumpImage(pCache->pDeflY,pCache->dimension,"deflY");
*/
	}

	/* change x and y into grid units in the cache which has same angular size as image pix */
	xgrid = x/pCache->pix_size+pCache->dimension[0]/2;
	ygrid = y/pCache->pix_size+pCache->dimension[1]/2;
	x1 = floor(xgrid);
	x2 = ceil(xgrid);
	if (x1==x2) {
		x2 +=1;
	}
	y1 = floor(ygrid);
	y2 = ceil(ygrid);
	if (y1==y2) {
		y2 +=1;
	}

	/* we now have the 4 points in the grid closest to the desired point */
	/* if these are outside the grid, then we have a small problem... */
	if(x1 < 0) {
		iDeltaX=-x1;
		x1=0;
	}
	if(x2 >= pCache->dimension[0]) {
		iDeltaX = (x2 - pCache->dimension[0]+1);
		x2=pCache->dimension[0]-1;
	}
	if(y1 < 0) {
		iDeltaY=-y1;
		y1=0;
	}
	if(y2 >= pCache->dimension[1]) {
		iDeltaY = (y2 - pCache->dimension[1]+1);
		y2=pCache->dimension[1]-1;
	}

	if (iDeltaX >0 || iDeltaY >0) {
		sprintf(strMessage,"====Warning: point exceeds range of cache. Desired: (%g,%g), size: (%d,%d)",xgrid,ygrid,(int)pCache->dimension[0],(int)pCache->dimension[1]);
		LOG_ERR(strMessage);
	}

	/* now interpolate between the closest points */
	t = (xgrid-x1)/(double)(x2-x1);
	u = (ygrid-y1)/(double)(y2-y1);
	*pDeflx = (1.0-t)*(1.0-u)*pCache->pDeflX[y1*pCache->dimension[0]+x1]
			+ t*(1.0-u)*pCache->pDeflX[y1*pCache->dimension[0]+x2]
			+ t*u*pCache->pDeflX[y2*pCache->dimension[0]+x2]
			+ (1.0-t)*u*pCache->pDeflX[y2*pCache->dimension[0]+x1];
	*pDefly = (1.0-t)*(1.0-u)*pCache->pDeflY[y1*pCache->dimension[0]+x1]
			+ t*(1.0-u)*pCache->pDeflY[y1*pCache->dimension[0]+x2]
			+ t*u*pCache->pDeflY[y2*pCache->dimension[0]+x2]
			+ (1.0-t)*u*pCache->pDeflY[y2*pCache->dimension[0]+x1];

EXIT:
	TRACE_OUT;
	return iRes;
}


/***************************
Function:	lm_CalcSersicDefl
Description: Entry point for calculation of a Sersic profile deflection
Arguments:
Returns:
****************************/
int	lm_CalcSersicDefl(real_t fx,real_t fy,real_t fAxRatio,real_t fScale,real_t fM, real_t *pDeflx,real_t *pDefly) {

	static	real_t	fLastM=-1;
	static	gsl_integration_workspace   *workspace=NULL;

	gsl_function    gfuncx,gfuncy;
	double	res=0, abserr=0;
    real_t params[1];
	int iResult = 0;

	TRACE_IN(lm_CalcSersicDefl);

	if (fx ==0 && fy==0) {
		*pDeflx=0;
		*pDefly=0;
		goto EXIT;
	}

	gfuncx.function = &lm_SersicDeflX_gsl;
	gfuncx.params = (void *)params;
	gfuncy.function = &lm_SersicDeflY_gsl;
	gfuncy.params = (void *)params;
	if (workspace==NULL) {
		workspace = gsl_integration_workspace_alloc(1000);
	}

	/* check for special case of circularly symmetric deVauc which
	 * has an analytical solution */
	if (fAxRatio == 1. && fM ==4.) {
		real_t	r,fCircDefl;

		r=sqrt(fx*fx + fy*fy);
		fCircDefl = lm_CalcDeVaucCircDefl(r,fScale);
		*pDeflx = fx/r*fCircDefl;
		*pDefly = fy/r*fCircDefl;
		goto EXIT;
	}

	if(fM != fLastM) {
		g_sersic_b = lm_calc_sersic_b(fM);
		fLastM = fM;
	}
    params[0] = 1./fM;

	g_axratio = fAxRatio;
	g_scale = fScale;
	g_tempx = fx;
	g_tempy = fy;

	iResult = gsl_integration_qag(&gfuncx,0.0,1.0,0.0,1e-6,1000,GSL_INTEG_GAUSS61,workspace,&res,&abserr);
    if (iResult != 0) {
        sprintf(strMessage,"ERROR: integration problems. return val: %d, x: %g, y:%g, axratio: %g, scale: %g, res: %g, abserr: %g\n",
            iResult,fx,fy,fAxRatio,fScale,res,abserr);
        LOG_ERR(strMessage);
		goto EXIT;
    }
	*pDeflx = fx*res;

	iResult = gsl_integration_qag(&gfuncy,0.0,1.0,0.0,1e-6,1000,GSL_INTEG_GAUSS61,workspace,&res,&abserr);
    if (iResult != 0) {
        sprintf(strMessage,"ERROR: integration problems. return val: %d, x: %g, y:%g, axratio: %g, scale: %g, res: %g, abserr: %g\n",
            iResult,fx,fy,fAxRatio,fScale,res,abserr);
        LOG_ERR(strMessage);
		goto EXIT;
    }
	*pDefly = fy*res;
/*
*	*pDeflx = fx*(qromo(lm_SersicDeflX,0,.01,midpnt)+qromo(lm_SersicDeflX,0.01,1,midpnt));
*	*pDefly = fy*(qromo(lm_SersicDeflY,0,.01,midpnt)+qromo(lm_SersicDeflY,0.01,1,midpnt));
*/

/*
printf("x,y: %g,%g. axratio: %g, scale: %g. defl: %g,%g.\n",fx,fy,fAxRatio,fScale,*pDeflx,*pDefly);
*/
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
/* old, numerical recipes version
*static real_t	lm_SersicDeflX(real_t x) {
*	double	t1,top;
*	t1 = 1.0 - (1.0 - g_axratio*g_axratio)*x;
*	top = x*(g_tempx*g_tempx + g_tempy*g_tempy/t1);
*	return exp(-g_sersic_b*pow(sqrt(top)/g_scale,))/sqrt(t1);
*}
*/
/* params: 1./M */
static double	lm_SersicDeflX_gsl(double x, void *params) {
	double	t1,top;
    real_t m_inv = *((real_t *)params);

	t1 = 1.0 - (1.0 - g_axratio*g_axratio)*x;
	top = x*(g_tempx*g_tempx + g_tempy*g_tempy/t1);
	return exp(-g_sersic_b*pow(sqrt(top)/g_scale,m_inv))/sqrt(t1);
}


/***************************
Function:
Description:
Arguments:
Returns:
****************************/
/* old numerical recipes version
*static	real_t	lm_SersicDeflY(real_t x) {
*	double	t1,top;
*	t1 = 1.0 - (1.0-g_axratio*g_axratio)*x;
*	top = x*(g_tempx*g_tempx + g_tempy*g_tempy/t1);
*	return exp(-g_sersic_b*pow(sqrt(top)/g_scale,0.25))/pow(t1,1.5);
*}
*/
static	double	lm_SersicDeflY_gsl(double x, void *params) {
	double	t1,top;
    real_t m_inv = *((real_t *)params);

	t1 = 1.0 - (1.0-g_axratio*g_axratio)*x;
	top = x*(g_tempx*g_tempx + g_tempy*g_tempy/t1);
	return exp(-g_sersic_b*pow(sqrt(top)/g_scale,m_inv))/pow(t1,1.5);
}


/***************************
Function:	lm_CalcDeVaucCircDefl
Description: Calculate de Vaucouleurs profile deflection for circularly symmetric lens
Arguments:	Re: DeVauc effective radius, r: position to calc deflection at
Returns:	unscaled (kappa=1) deflection
****************************/
static	real_t	lm_CalcDeVaucCircDefl(real_t r, real_t Re) {
	double	zeta;

	if (Re <= 0.) return 0;
	if (r ==0.) return 0;

	zeta = 7.66925*pow(r/Re,0.25);
	return 0.00336634*Re*(Re/r)*(1.0 - exp(-zeta)*(1+zeta*(1+zeta/2*(1+zeta/3*(1+zeta/4*(1+zeta/5*(1+zeta/6*(1+zeta/7))))))));
}


/***************************
Function:	lm_hernqf
Description: calculate the magic function for the Hernquist profile
Arguments:
Returns:
****************************/
static	real_t	lm_hernqf(real_t x) {
	real_t	temp=0.;

	if (x == 1.) {
		return 1.0;
	}
	if (x == 0.) {
		return 20.;
	}
	if (x < 1.) {
		temp = sqrt(1. - x*x);
		return lm_arctanh(temp)/temp;
	}
	if (x > 1.) {
		temp = sqrt(x*x -1.);
		return atan(temp)/temp;
	}
	return 1.0; 	/* doesn't get here, but shut the compiler up */
}

/***************************
Function:	lm_hernkappa
Description:	return the surface mass density (kappa) for a Hernquist
			profile without the central density included.
Arguments:	x: radial coordinate in units of scale length
Returns:	unscaled kappa
****************************/
real_t	lm_hernkappa(real_t x) {

	if (x > 0.99 && x < 1.01) {
		return 0.2665;
	}
	if (x < 1e-5) {
		return 20.0;
	}
	return (lm_hernqf(x)*(x*x + 2.) -3.)/((x*x-1.)*(x*x-1.));
}

/***************************
Function:	lm_CalcFerrersDefl
Description: calculate deflection angles for an l=2 ferrers ellipse.
			uses analytic expansions/recursion relations from Robert Schmidt's
			MSc thesis (University of Melbourne).
Arguments:	(x,y) position to calculate deflection for
			(a,b) major/minor axis of ellipse
Returns:	0: success, non-zero: failure
****************************/
int	lm_CalcFerrersDefl(real_t x, real_t y, real_t a, real_t b, real_t *fResX, real_t *fResY) {
	double	fRsqu=0, fLambda=0, a_squ, b_squ,x_squ, y_squ;

	*fResX = *fResY =0.0;

	if (a <= 0.0 || b <= 0.0) return 1;
	if (a == b) b = a*0.9999;
	a_squ = a*a;
	b_squ = b*b; 
	x_squ = x*x;
	y_squ = y*y;
	fRsqu = x_squ/a_squ + y_squ/b_squ;
	
	/* if we are outside the bar, then we must solve for lambda */
	if (fRsqu > 1.0) {
		double b1,c1;

		b1 = b_squ + a_squ - y_squ - x_squ;
		c1 = a_squ*b_squ - a_squ*y_squ - b_squ*x_squ;
		fLambda = 0.5*(-b1 + sqrt(b1*b1 - 4.0*c1));
	}

	/* calc the deflections */
	*fResX = a*b*x*(lm_ferrersq10(a_squ,b_squ,fLambda) + x_squ*(x_squ*lm_ferrersq30(a_squ,b_squ,fLambda) + 2.0*(y_squ*lm_ferrersq21(a_squ,b_squ,fLambda) - lm_ferrersq20(a_squ,b_squ,fLambda))) + y_squ*(y_squ*lm_ferrersq12(a_squ,b_squ,fLambda) - 2.0*lm_ferrersq11(a_squ,b_squ,fLambda)));
	*fResY = a*b*y*(lm_ferrersq01(a_squ,b_squ,fLambda) + x_squ*(x_squ*lm_ferrersq21(a_squ,b_squ,fLambda) + 2.0*(y_squ*lm_ferrersq12(a_squ,b_squ,fLambda) - lm_ferrersq11(a_squ,b_squ,fLambda))) + y_squ*(y_squ*lm_ferrersq03(a_squ,b_squ,fLambda) - 2.0*lm_ferrersq02(a_squ,b_squ,fLambda)));
	/*
	printf("defl for (a,b) (%g,%g) with (x,y) (%g,%g) is (%g,%g)\n",a,b,x,y,*fResX,*fResY);
	*/
	return 0;
}

/***************************
Function:
Description:
Arguments:
Returns:
****************************/
static	double	lm_ferrersq20(real_t a_squ, real_t b_squ, double l) {
	return (2.0/(sqrt(b_squ + l)*pow((a_squ+l),1.5)) - lm_ferrersq11(a_squ,b_squ,l))/3.0;
}

/***************************
Function:
Description:
Arguments:
Returns:
****************************/
static	double	lm_ferrersq02(real_t a_squ, real_t b_squ, double l) {
	return (2.0/(sqrt(a_squ + l)*pow((b_squ+l),1.5)) - lm_ferrersq11(a_squ,b_squ,l))/3.0;
}

/***************************
Function:
Description:
Arguments:
Returns:
****************************/
static	double	lm_ferrersq21(real_t a_squ, real_t b_squ, double l) {
	return (lm_ferrersq11(a_squ,b_squ,l) - lm_ferrersq20(a_squ,b_squ,l))/(a_squ - b_squ);
}

/***************************
Function:
Description:
Arguments:
Returns:
****************************/
static	double	lm_ferrersq12(real_t a_squ, real_t b_squ, double l) {
	return (lm_ferrersq02(a_squ,b_squ,l) - lm_ferrersq11(a_squ,b_squ,l))/(a_squ - b_squ);
}

/***************************
Function:
Description:
Arguments:
Returns:
****************************/
static	double	lm_ferrersq30(real_t a_squ, real_t b_squ, double l) {
	return (2.0/(sqrt(b_squ + l)*pow((a_squ+l),2.5)) - lm_ferrersq21(a_squ,b_squ,l))/5.0;
}

/***************************
Function:
Description:
Arguments:
Returns:
****************************/
static	double	lm_ferrersq03(real_t a_squ, real_t b_squ, double l) {
	return (2.0/(sqrt(a_squ + l)*pow((b_squ+l),2.5)) - lm_ferrersq12(a_squ,b_squ,l))/5.0;
}

/***************************
Function:
Description:
Arguments:
Returns:
****************************/
static	double	lm_ferrersq10(real_t a_squ, real_t b_squ, double lambda) {

	return 2.0*(1.0 - sqrt((b_squ + lambda)/(a_squ + lambda)))/(a_squ-b_squ);
}

/***************************
Function:
Description:
Arguments:
Returns:
****************************/
static	double	lm_ferrersq01(real_t a_squ, real_t b_squ, double lambda) {

	return 2*(sqrt((a_squ+lambda)/(b_squ+lambda)) -1.0)/(a_squ-b_squ);
}

/***************************
Function:
Description:
Arguments:
Returns:
****************************/
static	double	lm_ferrersq11(real_t a_squ, real_t b_squ, double lambda) {

	return (lm_ferrersq01(a_squ,b_squ,lambda) - lm_ferrersq10(a_squ,b_squ,lambda))/(a_squ-b_squ);
}

/***************************
Function:		lm_calc_sersic_b
Description:	Calculate the b factor for a Sersic light profile.
				From Ciotti & Bertin 1999
Arguments:		m: the Sersic function slope where I = exp(-b(r/Ro)^(1/m))
Returns:		the b factor
****************************/
static	real_t	lm_calc_sersic_b(real_t m) {
	return 2.0*m - (1.0/3.0) + 4.0/(405*m) + 46/(25515*m*m);
}

