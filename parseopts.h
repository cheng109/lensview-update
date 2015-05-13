#ifndef PARSEOPTS_H
#define PARSEOPTS_H
/* Copyright (C) 2006. Randall Wayth */

/**************************
Module:	parseopts.c
Author: Randall Wayth
Created:	4 Feb, 2000
Description: General purpose command line parsing utility
$Id: parseopts.h,v 1.5 2008/10/27 01:07:24 rwayth Exp rwayth $

Usage: 	This module provides a generic utility for parsing command line
		options. Options which have an associated value (e.g. "make -f filename")
		and options which activate a switch (e.g. "grep -i") are supported.
		Options which have an assocaited value can be an integer or a string.

		Command line options are parsed after examining the option string and
		variables are set accordingly. The function uses a variable argument list
		(like printf) and the interpretation of the variables depends on the
		contents of the option string.

		The option string (strOptionString) takes a format similar to printf
		except that option names are surrounded by delimiters (defined below)
		and followed by a "type" indicator of a single char. The general format
		of the string is:

		"<delimiter><option><delimiter><type><delimiter><option><delimiter><type>..."

		Option names can be more or less anything and don't have to have a "-"
		in front of them although this is good ettiquite. They cannot include
		the delimited character (%).

		Four types of options are supported:
		'i'-	integer (a number). Code just does an atoi() on this
		's'-	string. Just copies it as it.
		'f'-	floating point number. Just does an atof() on this
		't'-	toggle. Sets the associated variable to 1.
					( Does not check current value)

		Example option strings:
		"%-file%s"				- single option expecting a string
		"%-n%i"					- single option expecting a number
		"%-n%f"					- single option expecting a number with optional decimal point
		"%-quiet%t"				- single option which will set a flag to 1
		"%-quiet%t%-n%i%-file%s"	- combination of above
		"%-n%i%-quiet%t%-file%s"	- same again. Order doesn't matter

		All arguments in the variable arg list MUST be pointers. Strings must have
		space allocated already.

Example:
		Using the option string above.

		int	main(int argc, char *argv[]) {

			int	quiet_flag=0, num_widgets=0;
			char	filename[80]="";

			parse_options(argc, argv, "%-quiet%t%-n%i%-file%s", &quiet_flag, &num_widgets, filename);
			...
		}

**************************/

#define	PO_DELIMITER	'%'
#define	PO_MAX_ARG_LEN	80
#define	PO_INT			'i'
#define	PO_STRING		's'
#define	PO_TOGGLE		't'
#define	PO_FLOAT		'f'

/*************************
Function:	parse_options
Desctription:
Arguments:	iArgc:	argc from command line/main function
			pArgv:	argv[] from command line/main
			strOptionString:	formatted string which tells the
					function how to interpret the command line
			variable arg list:	variable argument list which the
					function will use depending on the contents of
					the option string
Returns:	0- SUCCESS
			1- FAILURE
**************************/

int	parse_options(int iArgc, char * const pArgv[], char *strOptionString, ...);

#endif /* ifndef PARSEOPTS_H */

