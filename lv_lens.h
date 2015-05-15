#ifndef	LV_LENS_H
#define LV_LENS_H
/******************
header files for the lv_lens module
of the "lensview" software
Copyright (C) 2006. Randall Wayth.
$Id: lv_lens.h,v 1.20 2008/10/27 01:07:24 rwayth Exp rwayth $
*******************/

/********************
Public defines
********************/
#define	LM_PTMASS			1
#define	LM_ISO_SPHERE		2
#define	LM_PIEP				3		/* elliptical lens potential with const density core */
#define	LM_SPEMD			4		/* softened elliptical power law model  */
#define	LM_NFW				5
#define	LM_SIE				6		/* elliptical surface mass density. Singular Isothermal ellipsoid */
#define	LM_EXTSHEAR			7
#define	LM_MASSSHEET		8
#define	LM_EXPDISC			9
#define	LM_SERSIC			10
#define	LM_FERRERS			11
#define	LM_USERDEF			12

#define	LM_LHMARK		"{"
#define	LM_RHMARK		"}"

#define	LM_NM_LENSCOMP	"lenscomp"
#define	LM_NAME_PTMASS	"PTMASS"
#define	LM_NAME_SIS		"SIS"
#define	LM_NAME_PIEP	"PIEP"
#define	LM_NAME_SPEMD	"SPEMD"
#define	LM_NAME_NFW		"NFW"
#define	LM_NAME_SIE		"SIE"
#define	LM_NAME_EXTSH	"EXTSHEAR"
#define	LM_NAME_MASSSHT	"MASSSHT"
#define	LM_NAME_EXPDISC	"EXPDISC"
#define	LM_NAME_SERSIC	"SERSIC"
#define	LM_NAME_FERRERS	"FERRERS"
#define LM_NAME_USERDEF	"USERDEF"

#define	LM_IGNORE_PROJ	42
#define LM_BADPROJ      -1

/********************
Public Structure definitions
********************/

/*******************
Public function definitions
*******************/
int			lm_CalcDeflection(lv_lensmodel_t *pLens, real_t fX, real_t fY, real_t *pDeltaX, real_t *pDeltaY, real_t *pMagnification);
lv_lensmodel_t	*lm_CreateLensModel(int	iNumCompoents);
void		lm_FreeLensModel(lv_lensmodel_t *pLensModel);
int			lm_InitLensModel(lv_lensmodel_t *pLens, int iType, real_t fXoffset, real_t fYoffset, int  iNumparams, real_t fParamsfrom[], real_t fParamsto[], real_t fParamsinc[], char *strNames[]);
int lm_CreateLMComp_ExtShear(lv_lensmodel_t *pLens, real_t fShearFrom, real_t fShearTo, real_t fShearInc,
    real_t fAnlgeFrom, real_t fAngleTo, real_t fAngleInc);
int lm_CreateLMComp_NFW(lv_lensmodel_t *pLens, real_t fXoffset, real_t fYoffset, real_t fMassScaleFrom,
	real_t fMassScaleTo, real_t fMassScaleInc, real_t fScaleLenFrom, real_t fScaleLenTo, real_t fScaleLenInc,
	real_t fEllipFrom, real_t fEllipTo, real_t fEllipInc, real_t fAngleFrom, real_t fAngleTo, real_t fAngleInc);
int lm_CreateLMComp_PtMass(lv_lensmodel_t *pLens, real_t fXoffset, real_t fYoffset, real_t fMassFrom, real_t fMassTo, real_t fMassInc);
int lm_CreateLMComp_MassSheet(lv_lensmodel_t *pLens, real_t fScaleFrom, real_t fScaleTo, real_t fScaleInc);
int lm_CreateLMComp_PIEP(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t fMassScaleFrom,
	real_t fMassScaleTo, real_t fMassScaleInc, real_t fEllipFrom, real_t fEllipTo, real_t fEllipInc,
	real_t fAngleFrom, real_t fAngleTo, real_t fAngleInc, real_t fCoreFrom, real_t fCoreTo, real_t fCoreInc);
int lm_CreateLMComp_SIE(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t fMassScaleFrom,
	real_t fMassScaleTo, real_t fMassScaleInc, real_t fEllipFrom, real_t fEllipTo, real_t fEllipInc,
			real_t fAngleFrom, real_t fAngleTo, real_t fAngleInc,
			real_t fCenterXFromm, real_t fCenterXTo, real_t fCenterXInc,
			real_t fCenterYFrom, real_t CenterYTo, real_t CenterYInc);
int lm_CreateLMComp_ExpDisc(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t fMassScaleFrom,
	real_t fMassScaleTo, real_t fMassScaleInc, real_t fEllipFrom, real_t fEllipTo, real_t fEllipInc,
	real_t fAngleFrom, real_t fAngleTo, real_t fAngleInc, real_t fScaleLenFrom, real_t fScaleLenTo, real_t fScaleLenInc);
int lm_CreateLMComp_Sersic(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t fMassScaleFrom,
	real_t fMassScaleTo, real_t fMassScaleInc, real_t fEllipFrom, real_t fEllipTo, real_t fEllipInc,
	real_t fAngleFrom, real_t fAngleTo, real_t fAngleInc, real_t fScaleLenFrom, real_t fScaleLenTo, real_t fScaleLenInc,
	real_t fMFrom, real_t fMto, real_t fMinc);
int lm_CreateLMComp_SIS(lv_lensmodel_t *pLens,real_t fXoffset, real_t fYoffset, real_t critfrom, real_t critto, real_t critinc);
int lm_CalcDeflComponent(lv_lenscomp *pLens, real_t fX, real_t fY, real_t *pDeltaX, real_t *pDeltaY);
lv_lensmodel_t *lm_ReadLensFile(char    *strFile);
int lm_decodeLensCompType(char *strName);
real_t	lm_expdiscdefly(real_t x, real_t y, real_t fAxRatio, real_t fScaleLen);
real_t	lm_expdiscdeflx(real_t x, real_t y, real_t fAxRatio, real_t fScaleLen);
real_t	lm_hernkappa(real_t x);
int lm_CalcSersicDefl(real_t fx,real_t fy,real_t fAxRatio,real_t fScale,real_t fM, real_t *pDeflx,real_t *pDefly);

#endif	/* ifndef LV_LENS_H */
