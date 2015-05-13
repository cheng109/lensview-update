/********************************
$Id: log.c,v 1.8 2006/11/20 14:37:04 rwayth Exp $
Generic Logging package for C programs

Copyright(C) Randall Wayth. Jan 2006.

Provides a generic logging tool which writes program messages to a log file, not stdout/stderr.
The log will be created in the directory pointed to by g_strLogFilePath (must be set).
A global variable defines the priority of messages to be written to the log. Messages with
a priority lower than this level are NOT written to the log.

To use this properly, your C functions should look like:
void myfunc() {
	int a,b,c;
	... other variables ...

	TRACE_IN(myfunc);

	... function code ...
	if (something bad happens) {
		printf("Something bad happened\n");
		LOG_ERR("Something bad happened");
		goto EXIT;
	}
	... function code ...

EXIT:
	TRACE_OUT;
}

There MUST be a TRACE_OUT call for all TRACE_INs. This forces the use of gotos if
the function must exit in the middle of code. Don't panic- it works fine.

Functions which do not have the TRACE_IN/TRACE_OUT combination should not use the TRACE() macro.
It will work, but the module name logged to the file will not be correct.
*********************************/

#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/utsname.h>
#include "log.h"
#include "common.h"

/*
 * Public Globals Vars
 */

int		g_iTracePriority = LOG_TRACE_NONE;
char	*g_strLogFilePath = NULL;
char	strMessage[MAX_MESG_LEN];			/* common string buffer for program to write messages into before tracing */

/*
 * Private Global Variables
 */

static	bool	bFileOpen = FALSE;
static	FILE	*traceFile=NULL;

/*
 * Private Function Prototypes
 */

static	void	log_OpenLogFile();

/*
 * Public Functions
 */
/***************************
Function:	log_LogTracePoint
Description:	Add a log entry to the log file if the incoming priority
				is <= the global level (in g_iTracePriority).
Arguments:	iPriority: trace priority (see log.h)
			strMessage: the message to log
			strModule: function/module name
Returns:	void
****************************/
void	log_LogTracePoint(int iPriority, char *strMessage,char *strModule) {

	char	strFullMessage[MAX_MESG_LEN + 100] = "";
	char	strDate[24] = "";
	int		iLen=0;
	time_t	now;

	/* file open yet? */
	if (bFileOpen == FALSE) {
		log_OpenLogFile();
	}

	/* only trace if this matches (or better) the priority */
	if (iPriority < g_iTracePriority) {
		return;
	}

	/* get the date/time */
	time(&now);
    strftime(strDate,sizeof(strDate),"%Y/%m/%d %H:%M:%S",localtime(&now));

	/* trim the message in case it's too long */
	if (strlen(strMessage) > (sizeof(strFullMessage) - sizeof(strDate) - strlen(strModule))) {
		strMessage[sizeof(strFullMessage)-1] = '\0';
	}

	/* generate the full line of the trace: date, module, message */
	sprintf(strFullMessage,"%s %s: %s",strDate, strModule, strMessage);

	/* make sure it has a \n at the end */
	iLen = strlen(strFullMessage);
	if ( (iLen > 0) && (strFullMessage[iLen-1] != '\n') )	{
		strFullMessage[iLen] = '\n';
		strFullMessage[iLen+1] = '\0';
	}

	fprintf(traceFile,strFullMessage);
	if (g_iTracePriority > LOG_LOW_PRI) {
		fflush(traceFile);
	}
}


/*
 * Private Functions
 */
/***************************
Function:	log_OpenLogFile
Description:	Create a log file. Uses the global variable g_strLogFilePath.
				Should be called when the first trace log is performed.
Arguments:	none
Returns:	void
****************************/
static	void	log_OpenLogFile() {

	char	*strFullFileName = NULL;
	char	strTime[50] = "";
	bool	bInsExtraSlash = FALSE;
	struct  utsname mach_name;
	time_t	now;

	/* set the file open to be true. Remains true even for failure */
	bFileOpen = TRUE;
	traceFile = NULL;

	/* get the machine name to include in the log file name */
	uname(&mach_name);

	/* create a log file with unique name based on the date when called */
	do {
		if (traceFile != NULL) {
			fclose(traceFile);
			fprintf(stderr,"WARNING: log file with name <%s> already exists. Waiting 1 second for a new name...\n",strFullFileName);
			sleep(1);
		}

		/* generate the log file name which includes the date */
		time(&now);
		strftime(strTime,sizeof(strTime),"%Y%m%d_%H%M%S",localtime(&now));

		if (strFullFileName == NULL) {
			strFullFileName = (char *) malloc(sizeof(char) * ((g_strLogFilePath == NULL? 0: strlen(g_strLogFilePath)) + strlen(strTime) + strlen(LOG_FILE_NAME_PREFIX) + strlen(mach_name.nodename) +11));
			if (strFullFileName == NULL) {
				fprintf(stderr,"log_OpenLogFile: no malloc\n");
				return;
			}
		}

		if (g_strLogFilePath != NULL && g_strLogFilePath[0] != '\0' && g_strLogFilePath[strlen(g_strLogFilePath)-1] != '/') {
			bInsExtraSlash = TRUE;
		}

		sprintf(strFullFileName,"%s%s%s_%s_%s.log",g_strLogFilePath == NULL? "": g_strLogFilePath,(bInsExtraSlash? "/": "" ), LOG_FILE_NAME_PREFIX,strTime,mach_name.nodename);
	}
	while ((traceFile = fopen(strFullFileName,"r")) != NULL);

	fprintf(stdout,"creating log file: %s\n",strFullFileName);

	traceFile = fopen(strFullFileName,"w+");
	if (traceFile == NULL) {
		fprintf(stderr,"WARNING: Could not open file name <%s>. No log will be generated.\n",strFullFileName);
		g_iTracePriority = LOG_TRACE_NONE;
	}

	if (strFullFileName != NULL) {
		free(strFullFileName);
	}
}
