#ifndef	LV_GEOM_H
#define	LV_GEOM_H

/********************
Header file for computational geometry functions
Copyright (C) 2006. Randall Wayth.
$Id: lv_geom.h,v 1.5 2006/09/17 21:58:51 rwayth Exp rwayth $
*********************/

#define	MAX_POLYGON_SIZE	10

/* placeholder for point types.*/
typedef	double	lv_pt_type_t;

typedef struct  _lv_point {
	lv_pt_type_t    x;
	lv_pt_type_t    y;
}	lv_point;

typedef	struct	_lv_line {
	lv_point	p1;
	lv_point	p2;
}	lv_line;

/* include extra space in the vertex array for sentinel values */
/* NOTE: Use vertex[0] and vertex[numVertices+1] for sentinel values */
/* do not use them for actual vertex points */
typedef	struct	_lv_polygon {
	int			numVertices;
	lv_point	vertex[MAX_POLYGON_SIZE+2];
}	lv_polygon;

/**********************
Public Funciton prototypes
**********************/

float	geomAreaOverlap(lv_polygon *p1, lv_polygon * p2);
int		geomIntersect(lv_line *pLine1, lv_line *pLine2);
int		geomIntersectPoint(lv_line *pLine1, lv_line *pLine2, lv_point *pResultPoint);
int		geomInsidePolygon(lv_point *pPoint, lv_polygon *pPolygon);
void	geomSortPolygon(lv_polygon *pPolygon);
float	geomPolygonArea(lv_polygon *pPolygon);

#endif	/* ifndef */
