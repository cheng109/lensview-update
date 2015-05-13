#ifndef	LV_CRITLINE_H
#define	LV_CRITLINE_H

/******************
 * header files for the lv_critline module
 * of the "lensview" software
 * Copyright (C) 2006. Randall Wayth.
 * $Id: lv_critline.h,v 1.3 2006/09/17 21:58:51 rwayth Exp rwayth $
 * *******************/

/********************
 * Public defines
 * ********************/

typedef struct _critline_pix_list {
		float	x1,x2;
		float	y1,y2;
		struct _critline_pix_list *pNext;
} cl_pix_list;

/* public function prototypes */
int cl_CalcCausticPoints(lv_lensmodel_t *pLens, cl_pix_list *pCritlineList, cl_pix_list *pCausticList);
int cl_CalcCritLinePoints(lv_lensmodel_t *pLens, lv_image_t  *pImg, cl_pix_list *pHead);
int cl_FreePixList(cl_pix_list *pList);
int cl_SavePixList(cl_pix_list *pList,char strFilename[]);
#endif /* ifndef LV_CRITLINE_H */
