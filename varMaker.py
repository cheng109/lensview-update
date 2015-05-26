"""
    @package
    @file numPhotons.py
    @count the tatal number of photons
    
    @brief Created by:
    @author Jun Cheng (Purdue)
    
    @warning This code is not fully validated
    and not ready for full release.  Please
    treat results with caution.
    
    Usage: python numPhotons.py image.fits num
    
    Notes: 1. 'num' means output 'peak +/- num' matrix
    
    """

import numpy as np
import sys
from astropy.io import fits

def readin(fileName):
   
    Xrange = matrix.shape[0]
    Yrange = matrix.shape[1]
    totalSum = np.sum(matrix)
    return totalSum, matrix, Xrange, Yrange

# Write out fits images for each order
def writeoutMatrix(inputFile, outputFile):
    
    
    hdulist = fits.open(inputFile)
    matrix = hdulist[0].data
    file = open(outputFile,'w')
    new_Matrix = matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if matrix[i][j]<50:
                new_Matrix[i][j]=50
            else: 
                new_Matrix[i][j] = matrix[i][j]
    output = fits.PrimaryHDU(new_Matrix)
    output.writeto(outputFile, clobber=True)

inputFile = sys.argv[1]
outputFile = sys.argv[2]

writeoutMatrix(inputFile,outputFile)

