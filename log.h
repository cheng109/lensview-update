#ifndef LOG_H
#define LOG_H
#include <limits.h>
/***************************
Generic Logging Package for C programs
Copyright (C) 2006. Randall Wayth.
$Id: log.h,v 1.6 2006/09/17 21:58:51 rwayth Exp rwayth $
***************************/

/*
 * Macro Definitions
 */

/* log priorities. */
#define LOG_TRACE_ALL	0
#define	LOG_LOW_PRI		1
#define	LOG_MED_PRI		2
#define	LOG_HIGH_PRI	3
#define	LOG_ERROR_PRI	4
#define	LOG_TRACE_NONE	5

#define	MAX_MESG_LEN	500
#define	LOG_FILE_NAME_PREFIX "lv"

/* macros actually used to create log entries */
#define	TRACE_IN(a)	{ char strLogFunctionName[] = #a; log_LogTracePoint(LOG_LOW_PRI,">",strLogFunctionName);
#define TRACE(a,b)	log_LogTracePoint((a),b,strLogFunctionName);
#define	TRACE_OUT	log_LogTracePoint(LOG_LOW_PRI,"<",strLogFunctionName); }
#define	LOG_ERR(a)	log_LogTracePoint(LOG_ERROR_PRI,a,strLogFunctionName);

/*
 * Function prototypes
 */

extern	void	log_LogTracePoint(int iPriority, char *strMessage, char *strModule);


/*
 * Global Variables
 */

extern	char	*g_strLogFilePath;
extern	int		g_iTracePriority;
extern	char	strMessage[MAX_MESG_LEN];
extern char    *g_paramFileName;
#endif /* LOG_H */
