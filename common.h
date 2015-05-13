/*****************
Copyright (C) 2006. Randall Wayth.
common definitions and macros for the lensview software
$Id: common.h,v 1.2 2006/09/17 21:58:51 rwayth Exp rwayth $
********************/
#ifndef	COMMON_H
#define	COMMON_H

/**************************
Common macro and type definitions
for C programs
**************************/

#ifndef TRUE
#define	TRUE	1
#endif

#ifndef FALSE
#define	FALSE	0
#endif

#define	MAX(a,b)	((a) < (b) ? (b) : (a))
#define	MIN(a,b)	((a) > (b) ? (b) : (a))

typedef int bool;

/***************************
Function:
Description:
Arguments:
Returns:
****************************/

#endif /* COMMON_H */
