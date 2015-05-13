/**************************
Functions for computational geometry.
Many functions based on Sedgewick, "Algorithms in C".
Copyright (C) Randall Wayth. 2006
$Id: lv_geom.c,v 1.9 2008/10/31 20:43:10 rwayth Exp rwayth $
*************************/

#include	<stdio.h>
#include	<stdlib.h>
#include	<math.h>
#include	<limits.h>
#include	<float.h>
#include	"log.h"
#include	"common.h"
#include	"lv_common.h"
#include	"lv_geom.h"

/***************************
Private function prototypes
****************************/
static	void	geomAddVertex(lv_point   *pPoint, lv_polygon *pPolygon);
static	int		geomCwDirection (lv_point   *pPoint0, lv_point *pPoint1, lv_point *pPoint2);
static  float   geomAngle(lv_point *pPoint1, lv_point *pPoint2);
static	float	geomCrossProd(lv_point   *pPoint0, lv_point *pPoint1, lv_point *pPoint2);

/***************************
Private Global Variables
***************************/


/***************************
Public funcitons
***************************/

/***************************
Funciton:	geomAreaOverlap
Description:	Compute the area overlap of a triangle (p2) of a square (p1). Will work
				for any two convex polygons
				The area of overlap can be calculated by finding the area of the
				polygon defined by a- the points of intersection of the polygons, and
				b- the points of each polygon which are inside the other
Args:		p1: pointer to polygon 1 (grid square)
			p2: pointer to polygon 2 (triangle)
Retruns:	area of overlap
****************************/
float	geomAreaOverlap(lv_polygon *p1, lv_polygon *p2) {
	bool	bInside;
	lv_line	l1,l2;
	lv_point	intersectPoint;
	lv_polygon	areaPolygon;
	int i,j;

	areaPolygon.numVertices =0;

	/* test the limiting cases, one is completely inside the other ?*/
	/* build a list of points of vertices inside the other polygons in
	 * case one isn't completely inside the other */

	/* is any of p1 inside p2? */
	bInside = TRUE;
	for (i=1; i <= p1->numVertices; i++) {
		j = geomInsidePolygon(&(p1->vertex[i]),p2);
		if (j==0) {
			bInside = FALSE;
		}
		else {
			/* add all of the points in p1 which are inside p2 */
			geomAddVertex(&p1->vertex[i],&areaPolygon);
		}
	}
	if (bInside) {
		/* p1 is completely inside p2 */
		geomSortPolygon(p1);
/*
printf("p1 completely inside p2. returning area of p1.\n");
*/
		return geomPolygonArea(p1);
	}

	/* is any of p2 inside p1? */
	bInside = TRUE;
	for (i=1; i<=p2->numVertices; i++) {
		if (geomInsidePolygon(&p2->vertex[i],p1)==0) {
			bInside = FALSE;
		}
		else {
			geomAddVertex(&p2->vertex[i],&areaPolygon);
/*
printf("Inside p1. point %f,%f\n",p2->vertex[i].x,p2->vertex[i].y);
*/
		}
	}
	if (bInside) {
		/* p2 is completely inside p1 */
		geomSortPolygon(p2);
/*
printf("p2 completely inside p1. returning area of p2.\n");
*/
		return geomPolygonArea(p2);
	}

	/* build a list of all the intersections of the edges of the polygons */
	/* geomAddVertex ignores ones we have already taken into account */
	for (i=1; i<= p1->numVertices; i++) {
		l1.p1 = p1->vertex[i];
		l1.p2 = p1->vertex[(i%p1->numVertices) + 1];
		for (j=1; j<= p2->numVertices; j++) {
			l2.p1 = p2->vertex[j];
			l2.p2 = p2->vertex[(j%p2->numVertices) +1];
			if (geomIntersectPoint(&l1,&l2,&intersectPoint)) {
				geomAddVertex(&intersectPoint,&areaPolygon);
			}
		}
	}

	geomSortPolygon(&areaPolygon);

/*
printf("Area overlap polygon.\n");
print_polygon(&areaPolygon);
*/

	return geomPolygonArea(&areaPolygon);
}

/***************************
Function:	geomPolygonArea
Description: Compute the area of a convex polygon. Assumes the polygon has been
			previously "sorted" so that the vertices are consecutive
			around the perimeter of the polygon.
Arguments:	pPolygon: pointer to the polygon
Returns:	the area
****************************/
float   geomPolygonArea(lv_polygon *pPolygon) {
	int		i;
	float	fArea=0;
	lv_point	*base;

	base = &(pPolygon->vertex[1]);

	/* divide the polygon into trangles and add up their area
		by taking the cross product of the edge vectors */
	for (i=3; i<=pPolygon->numVertices; i++) {
		fArea += 0.5 * fabs(geomCrossProd(base, &(pPolygon->vertex[i-1]), &(pPolygon->vertex[i])));
	}
	return fArea;
}

/***************************
Function:	geomSortPolygon
Description: sort the vertices of the polygon so that they form a
			consecutive series around the perimeter of the polygon
			with the first vertex (vertex[1], not 0) being least x/y val
Arguments:	pPolygon: pointer to the polygon
Returns:	none
****************************/
void	geomSortPolygon(lv_polygon *pPolygon) {

	int		i=0,j=0,iMinpos=0;
	lv_pt_type_t	yMin = FLT_MAX, xMin = FLT_MAX;
	lv_point	*base, temp;
	float		fMinAngle=FLT_MAX, fCurrent=0;

	/* find the smallest x-coord for the smallest y-coord. This is used
		by the geomInsidePolygon function as a starting point */
	for (i = 1; i<= pPolygon->numVertices; i++) {
		if (pPolygon->vertex[i].y <= yMin) {
			if (pPolygon->vertex[i].y < yMin) {
				/* if this is the smallest y so far, reset Xmin and remember yMin */
				yMin = pPolygon->vertex[i].y;
				xMin = FLT_MAX;
			}
			if (pPolygon->vertex[i].x < xMin) {
				/* if this is equal to the yMin, then check the xMin value */
				xMin = pPolygon->vertex[i].x;
				j=i;
			}
		}
	}

	/* swap the value of the smallest point into vertex[1] if it's not already there */
	if (j != 1) {
		temp = pPolygon->vertex[1];
		pPolygon->vertex[1] = pPolygon->vertex[j];
		pPolygon->vertex[j] = temp;
	}
	base = &(pPolygon->vertex[1]);

	/* selection sort the remaining vertices according to angle from the 1st vertex */
	/* sort in ascending order */
	for (i=2; i < pPolygon->numVertices; i++) {
		fMinAngle = geomAngle(base, &pPolygon->vertex[i]);
		iMinpos = i;
		for (j=i+1; j <= pPolygon->numVertices; j++) {
			fCurrent = geomAngle(base,&pPolygon->vertex[j]);
			/* remember the smallest angle of the current and the min */
			if (fCurrent < fMinAngle) {
				fMinAngle = fCurrent;
				iMinpos = j;
			}
		}

		/* swap the smallest into place */
		if (iMinpos != i) {
			temp = pPolygon->vertex[i];
			pPolygon->vertex[i] = pPolygon->vertex[iMinpos];
			pPolygon->vertex[iMinpos] = temp;
		}
	}
	/* lastly, set the sentinel values at position 0 and N+1 */
	pPolygon->vertex[0] = pPolygon->vertex[pPolygon->numVertices];
	pPolygon->vertex[pPolygon->numVertices+1] = pPolygon->vertex[1];
}

/***************************
Function:	geomInsidePolygon (from Sedgewick p354)
Description: Determine if a point is inside a polygon
			NOTE: Assumes points in polygon are sorted with vertex[1]
			is the smallest x co-ord among points with the smallest y coord
			see Sedgewick for details.
Arguments:	pPoint: The point. pPolygon: the polygon
Returns:	1- Yes, 0- no. Returns 1 for a point lying on the edge of the polygon
****************************/
int geomInsidePolygon(lv_point *pPoint, lv_polygon *pPolygon) {
	int	i=0, iCount=0, j=0, iDir =0, iCurrDir=0;
	lv_line	lt,lp;

	pPolygon->vertex[0] = pPolygon->vertex[pPolygon->numVertices];
	pPolygon->vertex[pPolygon->numVertices+1] = pPolygon->vertex[1];

/*	special case: for 3/4 sided convex polygon, we can use the geomCwDirection
 *	funciton to determine the result */
/* for a point inside a polygon, the direction between successive vertices should be
 * the same. For outside, there will be at least one direction reversal. For points on
 * the line, there is no direction, so ignore that result */
	if (pPolygon->numVertices > 2 && pPolygon->numVertices < 5) {
/*
printf("using special 3/4 sided convex polygon code.\n");
*/
		for (i=1; i <= pPolygon->numVertices; i++) {
			iCurrDir = geomCwDirection(pPoint, &pPolygon->vertex[i],&pPolygon->vertex[i+1]);
			if (iDir != 0 && iDir * iCurrDir == -1) {
				return 0;
			}
			else {
				if (iDir==0) {
					iDir = iCurrDir;
				}
			}
		}
		return 1;
	}

	lt.p1 = *pPoint; lt.p2 = *pPoint; lt.p2.x = INT_MAX;
	for (i=1; i<=pPolygon->numVertices; i++) {
		lp.p1 = pPolygon->vertex[i];
		lp.p2 = pPolygon->vertex[i];
		if (!geomIntersect(&lp,&lt)) {
			lp.p2 = pPolygon->vertex[j];
			j=i;
			if (geomIntersect(&lp,&lt)) {
				iCount++;
			}
		}
	}
	return iCount & 1;
}

/***************************
Function:	geomIntersectPoint
Description:	Test for intersection of two lines and return the point of intersection.
Arguments:
Returns:	0:	no intersection
			1:	intersection, pResultPoint updated
****************************/
int	geomIntersectPoint(lv_line *pLine1, lv_line *pLine2, lv_point *pResultPoint) {

	float	dx1, dx2, dy1, dy2, m2;

	if (!geomIntersect(pLine1,pLine2)) {
		return 0;
	}
	dx1 = pLine1->p1.x - pLine1->p2.x;
	dy1 = pLine1->p1.y - pLine1->p2.y;

	dx2 = pLine2->p1.x - pLine2->p2.x;
	dy2 = pLine2->p1.y - pLine2->p2.y;

	/* special case: line 1 is vertical */
	if (dx1 == 0) {
		if (dx2 == 0) {
			return 0;
		}
		else {
			pResultPoint->x = pLine1->p1.x;
			m2 = dy2/dx2;
			pResultPoint->y = m2* (pLine1->p1.x - pLine2->p1.x) + pLine2->p1.y ;
			return 1;
		}
	}

	/* special case: line 1 is horizontal  */
	if (dy1 == 0) {
		if (dy2 == 0) {
			return 0;
		}
		else {
			pResultPoint->y = pLine1->p1.y;
			m2 = dx2/dy2;
			pResultPoint->x = m2* (pLine1->p1.y - pLine2->p1.y) + pLine2->p1.x ;
			return 1;
		}
	}
	return 0; /* never reach here, but keep compiler happy */
}

/***************************
Function:	geomIntersect (from Sedgewick p 351)
Description: find if two lines intersect
Arguments:	pLine1, pLine2: pointers to line stucts
Returns:	1: yes
			0: no
****************************/
int geomIntersect(lv_line *pLine1, lv_line *pLine2) {

	/* check for special case: the two endpoints on a line are actually the same point. */
	if (pLine1->p1.x == pLine1->p2.x && pLine1->p1.y == pLine1->p2.y) {
	}

	if (pLine2->p1.x == pLine2->p2.x && pLine2->p1.y == pLine2->p2.y) {
	}

	return (( geomCwDirection(&pLine1->p1, &pLine1->p2, &pLine2->p1) *
			geomCwDirection(&pLine1->p1, &pLine1->p2, &pLine2->p2)) <= 0) &&
			(( geomCwDirection(&pLine2->p1, &pLine2->p2, &pLine1->p1) *
			geomCwDirection(&pLine2->p1, &pLine2->p2, &pLine1->p2)) <= 0);
}


/***************************
Private functions
***************************/

/***************************
Function:	geomAddVertex
Description:	add a vertex to a generic polygon. Ignores point if it already
				exists in the polygon
Arguments:	pPoint:	vertex point to be added
			pPolygon: polygon to add point to
Returns:	none
****************************/
static	void	geomAddVertex(lv_point   *pPoint, lv_polygon *pPolygon) {
	int		k;
	bool	bPointIncluded;

	/* check to see if this point is already included. This is possible
	 * if intersections are in vertices. Allow for small numerical differences */
	bPointIncluded = FALSE;
	for (k=1; k <= pPolygon->numVertices; k++) {
		if (fabs((pPolygon->vertex[k].x - pPoint->x)/(pPoint->x==0? 1e-8: pPoint->x)) < 1e-7 &&
					fabs((pPolygon->vertex[k].y - pPoint->y)/(pPoint->y==0? 1e-8: pPoint->y)) < 1e-7) {
			bPointIncluded = TRUE;
			break;
		}
	}
	if (bPointIncluded == FALSE) {
		pPolygon->numVertices++;
		pPolygon->vertex[pPolygon->numVertices] = *pPoint;
	}
}

/***************************
Function:	geomCwDirection (from Sedgewick page 350)
Description:	Determine whether traversing a set of points is in a
				clockwise or anticlockwise direction
Arguments:	pPoint0, pPoint1, pPoint2:	pointers to the points
Returns:	1:	anticlockwise
			-1:	clockwise
			0:	point 2 is on the line segment between point 0 and point 1
****************************/
static	int	geomCwDirection (lv_point	*pPoint0, lv_point *pPoint1, lv_point *pPoint2) {
	lv_pt_type_t	dx1, dx2, dy1, dy2, temp;

	dx1 = pPoint1->x - pPoint0->x;
	dy1 = pPoint1->y - pPoint0->y;
	dx2 = pPoint2->x - pPoint0->x;
	dy2 = pPoint2->y - pPoint0->y;

	temp = dy1*dx2 - dx1*dy2;
	if (temp < 0) return 1;
	if (temp > 0) return -1;
	if (dx1 == 0 && dx2 == 0) return 0;
	if (dy1 == 0 && dy2 == 0) return 0;
	if ( (dx1 * dx2 < 0) || (dy1 * dy2 < 0) ) return -1;
	if ( (dx1*dx1 + dy1*dy1) < (dx2*dx2 + dy2*dy2) ) return 1;
	return 0;
}

/***************************
Function:	geomAngle (from Sedgewick p 353)
Description: determine the relative angle between two points. NOTE: This function uses
			an approximation and doesn't return the "true" angle using arctan, but it can be used
			to determine if	one angle is greater than the other. Will be useful for sorting
			polygon vertices
Arguments:	pPoint1, pPoint2: points to find the angle between
Returns:	angle
****************************/
static	float	geomAngle(lv_point *pPoint1, lv_point *pPoint2) {
	lv_pt_type_t	dx, dy, ax, ay;
	float	fTangent=0;

	dx = pPoint2->x - pPoint1->x;
	ax = fabs(dx);
	dy = pPoint2->y - pPoint1->y;
	ay = fabs(dy);
	fTangent = (ax+ay == 0) ? 0 : (float) dy/(ax+ay);
	if (dx < 0) {
		fTangent = 2-fTangent;
	}
	else if (dy < 0) {
		fTangent += 4;
	}
	return fTangent * 90;
}

/***************************
Function:	geomCrossProd
Description: return the value of the cross product of two lines
			assumed to be in the X-Y plane. Point 0 is the common point
Arguments:	pPoint0, pPoint1, pPoint2: pointers to points.
Returns:	value of the cross product. Can be negative
****************************/
static  float   geomCrossProd(lv_point   *pPoint0, lv_point *pPoint1, lv_point *pPoint2) {
	return (pPoint1->x - pPoint0->x) * (pPoint2->y - pPoint0->y) -
			(pPoint2->x - pPoint0->x) * (pPoint1->y - pPoint0->y); 
}

