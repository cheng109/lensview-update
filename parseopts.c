/**************************
Module:	parseopts.c
Author:	Randall Wayth
Created:	4 Feb, 2000
Description: See parseopts.h for details
Copyright (C) 2006. Randall Wayth
$Id: parseopts.c,v 1.4 2008/10/14 15:11:54 rwayth Exp rwayth $
**************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "parseopts.h"

int	parse_options(int iArgc, char * const pArgv[], char *strOptionString, ...) {

	va_list	pVarArg;
	int		iArgIndex = 0, iResult = 0, *pVaNumArg=NULL, iNumDelim=0, i;
	float	*pVaFloatArg = NULL;
	char	*pOption = NULL, *pVaStrArg=NULL, *pScan = NULL;
	char	*pArg = NULL, strSearchArg[PO_MAX_ARG_LEN] = "";

	/* sanity check */
	if (strOptionString == NULL || pArgv == NULL) {
		fprintf(stderr,"parse_options: option list is NULL\n");
		return 1;
	}

	/* iterate through the command line and set flags/variables accordingly */
	/* since there is no restriction on the order of args, then we need
		to scan the whole option array each time. This isn't efficient, but
		it's a one-off so no big deal and we need to make sure all command
		line args are examined anyway. */
	iArgIndex = 1;				/* argv[0] is the name of the program */
								/* argc is 1 for a program with no args */

	while (iArgIndex < iArgc && iResult ==0) {

		/* init the variable arg list pointer. The last named arg must be used here. */
		va_start(pVarArg, strOptionString);

		pArg = pArgv[iArgIndex];
		pOption = NULL;
		/* scan the option string for this arg with the delimiter */
		/* if it's not there, then fail */

		sprintf(strSearchArg,"%c%s%c",PO_DELIMITER, pArg, PO_DELIMITER);

		pOption = strstr(strOptionString, strSearchArg);
		if (pOption == NULL) {
			fprintf(stderr,"Unknown option: %s\n",pArg);
			iResult = 1;
		}
		else {
			/* find the type of option */
			/* point to the type indicator for this option */
			pOption += strlen(strSearchArg);

			/* find how many options along the option list we are and
				move the same number of args through the variable arg
				list */
			/* move through the string up to where pOption is pointing.
				 the number of args to skip will be
				(the number of delimiters/2)-1 */
			iNumDelim =0;
			for (pScan = strOptionString; pScan != pOption; pScan++) {
				if (*pScan==PO_DELIMITER) {
					iNumDelim++;
				}
			}

			/* skip the required number of args. It doesn't matter if we get the
				type wrong here as long as all pointers are the same size */
			for (i=0; i<(iNumDelim/2)-1; i++) {
				pVaNumArg = va_arg(pVarArg, int *);
			}

			switch(*pOption) {
				case PO_TOGGLE: 
							pVaNumArg = va_arg(pVarArg, int *);
							*pVaNumArg = 1;
							iArgIndex++;
							break;
				case PO_STRING:
							pVaStrArg = va_arg(pVarArg, char *);
							strcpy(pVaStrArg,pArgv[iArgIndex+1]);
							iArgIndex += 2;
							break;
				case PO_INT:
							pVaNumArg = va_arg(pVarArg, int *);
							*pVaNumArg = atoi(pArgv[iArgIndex+1]);
							iArgIndex += 2;
							break;
				case PO_FLOAT:
							pVaFloatArg = va_arg(pVarArg, float *);
							*pVaFloatArg = (float) atof(pArgv[iArgIndex+1]);
							iArgIndex += 2;
							break;
				default:
							fprintf(stderr,"Unknown option type <%c>\n",*pOption);
							iResult =1;
							break;
			}
		}
		/* after processing, must do this */
		va_end(pVarArg);
	}

	return iResult;
}

