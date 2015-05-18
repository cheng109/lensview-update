import sys, os
import re

class Model():
    name = ""
    xoffset = 0;
    yoffset = 0;

    def _init_(self):
        this->name=name


def createWorkDirectory(dirName):
    # create a work directory:
    currentDir = os.getcwd()
    if not os.path.exists(currentDir+ "/"+ dirName):
        print "hello"
        os.mkdir("work")
    return


def readCommandFile(commandFile):

    infile = open(commandFile, 'r')
    for line in infile.readlines():
        m = re.match('.*paramfile\s+(\S+)\s+',line)

        if m:
            paramFileName = m.group(1);
            return paramFileName

    print "Parameter file name not found!"


def createMultipleParamFiles(paramFileName):
    infile = open(paramFileName, 'r')
    for line in infile.readlines():
        if line[0]!="#":
            temp = re.split(',|\{|\}|\(|\)|\n ',line)

            #temp = line.split(",")
            print temp
    return



def main():

    createWorkDirectory("work");
    paramFileName = readCommandFile("run")
    createMultipleParamFiles("example/"+paramFileName)


if __name__=='__main__':
    main()
    
