UPDATEs from Jun Cheng
======================
Models modified:  
- SIE
- NFW
- SPEMD
- PTMASS
- PIEP

New fetures: 
- Output the residual maps
- Output model and source images with names associated with parameter files ( Ready for Parallel computing ) 
- 'jobDistribution.py' can be used to split the jobs into small pieces and submit the PBS script to computing clusters



Lensview: Software for modelling resolved gravitational lenses.
===============================================================
Copyright (C) 2006-2008. Randall Wayth.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

PREREQUISITES
=============

Lensview depends on:
The GNU scientific library (GSL)
The CFITSIO library
The FFTW library
(optionally) the fastell library

C and F77 compilers.

See http://cfa-www.cfa.harvard.edu/~rwayth/lensview
for details.

INSTALLATION
============

- create a directory for lensview
- untar the tarball into this directory
- install the dependencies
- edit the makefile so that includes and libraries can be found by the compiler
- make lensview

Acknowledgments
===============
Lensview makes use of Rennan Barkarna's "fastell" code.
http://wise-obs.tau.ac.il/~barkana/fastell.f
fastell makes use of some of the slatec packages available from netlib:
http://www.netlib.org/slatec/

