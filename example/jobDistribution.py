import sys, os
import re
import itertools
import subprocess
import shutil
class Model:
    def __init__(self, line):
        self.fromList = []
        self.toList = []
        self.incList = []
        self.name=line[0]
        self.xoffset= float(line[1])
        self.yoffset= float(line[2])
        for i in range(3, len(line),3):
            self.fromList.append(float(line[i]))
        for i in range(4, len(line),3):
            self.toList.append(float(line[i]))
        for i in range(5, len(line),3):
            self.incList.append(float(line[i]))


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
        m = re.match('(.*paramfile)\s+(\S+)\s+(.*)',line)

        if m:
            paramFileName = m.group(2)

            return paramFileName, (m.group(1), m.group(3))

    print "Parameter file name not found!"

def sperateRegions(model, numFiles, upLimit):
    fromList = model.fromList
    toList = model.toList
    incList = model.incList
    NestList = []
    for i in range(len(fromList)):
        newInc = max((toList[i]-fromList[i])/numFiles, incList[i])
        start = fromList[i]
        temp = []
        while start<=toList[i]:
            if i<upLimit:
                temp.append((start, start+newInc, incList[i]))
                start+=newInc
            else:
                temp.append((start, toList[i], incList[i]))
                break


        NestList.append(temp)
    cross = list(itertools.product(*NestList))


    newParamFileList = [""]*len(cross)
    for i in range(len(cross)):
        tempString = "lenscomp{" + model.name + "}("

        for (a,b,c) in cross[i]:
            tempString += (str(round(a,2)) + "," + str(round(b,2)) + "," + str(round(c,2)) + ",")
        tempList = list(tempString)
        tempList[-1] = ")"
        newParamFileList[i] = "".join(tempList)+"\n"
    return newParamFileList


def writeOutParametersFile(cross, paramFileName):
    paramFileNameList = []
    for i in range(len(cross)):
        newFileName = paramFileName + "_" + str(1000000+i)
        paramFileNameList.append(newFileName)
        f = open(newFileName, 'w')
        for string in cross[i]:
            f.write(string)
        f.close()
    return paramFileNameList

def createScriptPBS(remain, paramFileNameList, commandFileName):
    head = """
hello! \n
"""
    for paramName in paramFileNameList:
        f = open(commandFileName+"_"+ paramName, 'w')
        f.write(head + remain[0] + " " +  paramName + " " + remain[1] + "\n" )
        f.close()
    return


def createMultipleParamFiles(paramFileName, numFiles, upLimit):
    modelList = []
    infile = open(paramFileName, 'r')
    for line in infile.readlines():
        if line[0]!="#":
            temp = re.split(',|\{|\}|\(|\)|\n',line)
            temp = [x for x in temp if (x!='' and x!='lenscomp')]
            m = Model(temp)
            modelList.append(m)
    piece = []
    for model in modelList:
        piece.append(sperateRegions(model, numFiles, upLimit))
    # Combine parameters from those models:
    cross = list(itertools.product(*piece))
    paramFileNameList = writeOutParametersFile(cross, paramFileName)
    return paramFileNameList




def main():
    sep = 2
    upLimit = 4
    commandFileName = sys.argv[1]
    workDirectoryName = "work"
    createWorkDirectory(workDirectoryName)
    shutil.copy(commandFileName, workDirectoryName+"/")
    paramFileName, remain = readCommandFile(commandFileName)
    shutil.copy(paramFileName, workDirectoryName+"/")
  # copy 'command' and 'parameter' files to work directory
    os.chdir(workDirectoryName)
    paramFileNameList = createMultipleParamFiles(paramFileName, sep, upLimit)
    createScriptPBS(remain, paramFileNameList, commandFileName)

if __name__=='__main__':
    main()
    
