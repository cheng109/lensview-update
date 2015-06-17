#!/usr/local/Ureka/variants/common/bin/python

import glob, os, sys
import pyfits
import operator

def getResidual(fileList, varFile):
    normRes = {}
    for fits in fileList:
        hdulist = pyfits.open(fits)
        data = hdulist[0].data
        sum = 0
        for i in range(data.shape[0]):
            for j in range(data.shape[1]):
                sum += abs(data[i][j])
        normRes[fits] = sum

    #normRes = sorted(normRes.items(), key=operator.itemgetter(1))
    return normRes


def writeOutRes(outputName, normRes):
    f = open(outputName,'w')
    f.write("#0 File name \n"+
            "#1 Sum of abs(residual)"+
            "\n\n")

    for key, value in normRes.iteritems():
        f.write(key + "\t" + str(value) + "\n")
    f.write("\nThe one with Minimum abs(residual) is:\n"+ str(min(normRes, key=normRes.get)) + " ")
    f.write(str(normRes[min(normRes, key=normRes.get)]) + "\n")
    f.close()



def main():
    workDir = sys.argv[1]
    #currentDir = "/Users/cheng109/work/LensView_update/example/"
    currentDir = os.path.dirname(os.path.abspath(__file__))
    print currentDir+ workDir
    resFileList=glob.glob(currentDir+ "/"+ workDir + "model_res*.fits")
    varFile=glob.glob(currentDir+ "/"+ workDir + "var4.fits")
    normRes = getResidual(resFileList, varFile)
    writeOutRes(workDir+ "output.txt", normRes)

if __name__=='__main__':
    main()
