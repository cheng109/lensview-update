/*************************

C headers for fastell.f fortran code.
Randall Wayth. 31 May, 2001.

To compile: Create object code files:

prompt> f90 -c -fast fastell.f dqnc79.f
prompt> cc -c <cprograms.c>

then link with C or Fortran compiler. Using the Fortran one is easier
since you don't have to worry about which libraries fortran uses.

prompt> cc -o myprog fastell.o cprogs.o dqnc79.o -lF77 -lm -lM77 -lsunmath		# Solaris
OR
prompt>	cc -o myprog fastell.o cprogs.o dqnc79.o -lm -lfor				# Digital OSF/1
OR
prompt> f90 -o myprog fastell.o cprogs.o dqnc79.o					# Solaris
OR
prompt> f90 -nofor_main -o myprog fastell.o cprogs.o dqnc79.o				# Digital OSF/1

If you're having trouble, the fortran man page will list the libraries
which it needs.

For more help on using fortran in C code, go to:
http://www.ictp.trieste.it/~manuals/programming/sun/fortran/prog_guide/11_cfort.doc.html

****************************/
/*
 *	fastell.f
 */

/* calculate the deflection angle using straight integration
 * and the fast way */
extern  void ellipdefl_(double *x1, double *x2, double *q, double *gamma, double *axisratio, double *coreradsqu, double deflection[2]);
extern  void fastelldefl_(double *x1, double *x2, double *q, double *gamma, double *axisratio, double *coreradsqu, double deflection[2]);

/* calculate the magnification and deflection using straight integration
 * and the fast way
 */
extern  void ellipmag_(double *x1, double *x2, double *q, double *gamma, double *axisratio, double *coreradsqu, double deflection[2],double mag[2][2]);
extern  void fastellmag_(double *x1, double *x2, double *q, double *gamma, double *axisratio, double *coreradsqu, double deflection[2],double mag[2][2]);

/* calculate the lens potential phi */
extern  void ellipphi_(double *x1, double *x2, double *q, double *gamma, double *axisratio, double *coreradsqu, double *phi);
