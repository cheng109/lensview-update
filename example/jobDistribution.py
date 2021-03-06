#!/usr/bin/python
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
        os.mkdir(dirName)
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
    quota = upLimit
    for i in range(len(fromList)):
        newInc = max((toList[i]-fromList[i])/numFiles, incList[i])
        start = fromList[i]
        temp = []
        if i<=quota and (fromList[i]==toList[i] or fromList[i]+newInc > toList[i]):
            quota +=1
            temp.append((fromList[i], toList[i], incList[i]))
        elif i>quota:
            temp.append((fromList[i], toList[i], incList[i]))
        else:
            while start+newInc<=toList[i]:
                temp.append((start, start+newInc, incList[i]))
                start+=newInc
            if start<toList[i] and i<upLimit:
                temp.append((start, toList[i], incList[i]))
        NestList.append(temp)
    cross = list(itertools.product(*NestList))

    newParamFileList = [""]*len(cross)
    for i in range(len(cross)):
        tempString = "lenscomp{" + model.name + "}(0,0," 
        for (a,b,c) in cross[i]:
            tempString += (str(round(a,2)) + "," + str(round(b,2)) + "," + str(round(c,2)) + ",")
        tempList = list(tempString)
        tempList[-1] = ")"
        newParamFileList[i] = "".join(tempList)+"\n"
    return newParamFileList


def writeOutParametersFile(cross, paramFileName):
    paramFileNameList = []
    for i in range(len(cross)):
        newFileName = paramFileName + "_ID_" + str(1000000+i)
        paramFileNameList.append(newFileName)
        f = open(newFileName, 'w')
        for string in cross[i]:
            f.write(string)
        f.close()
    return paramFileNameList

def createScriptPBS(remain, paramFileNameList, commandFileName):
    fullWorkPath = os.path.dirname(os.path.realpath(__file__))+"\n"
    head = """#PBS -q standby
#PBS -l nodes=1:ppn=1,naccesspolicy=shared
#PBS -l walltime=4:00:00
cd """
    head += fullWorkPath
    pbsFileNameList = []
    for paramName in paramFileNameList:
        fileName = commandFileName+"_"+ paramName
        pbsFileNameList.append(fileName)
        f = open(fileName, 'w')
        f.write(head + remain[0] + " " + paramName + " " + remain[1] + "\n" )
        f.close()
    return pbsFileNameList


def createMultipleParamFiles(paramFileName, numFiles, upLimit):
    modelList = []
    infile = open(paramFileName, 'r')
    for line in infile.readlines():
        if line[0]!="#":
            temp = re.split(',|\{|\}|\(|\)|\n',line)
            temp = [x for x in temp if (x!='' and x!='lenscomp')]
            if len(temp)>0:
                m = Model(temp)
                modelList.append(m)
    piece = []
    for model in modelList:
        piece.append(sperateRegions(model, numFiles, upLimit))
    # Combine parameters from those models:
    cross = list(itertools.product(*piece))
    paramFileNameList = writeOutParametersFile(cross, paramFileName)
    return paramFileNameList

def copyFiles(command, destination):
    for files in os.listdir("."):
        if files.endswith(".fits") or files.endswith(".fits.gz"):
            shutil.copy(files, destination)
    shutil.copy(command, destination)
    shutil.copy("lensview", destination)

def submitJobs(pbsFileNameList, machine):
    for fileName in pbsFileNameList:
        if machine=='server':
            subprocess.call(["qsub", fileName])
        elif machine=='local':
            subprocess.call(["bash", fileName])
        else:
            print "Please select 'server' or 'local'"
    return 0


def main():
    sep = 3
    upLimit = 3
    commandFileName = sys.argv[1]
    if len(sys.argv)==3:
        workDirectoryName = sys.argv[2]
    else: 
        workDirectoryName = "work"
    createWorkDirectory(workDirectoryName)
    copyFiles(commandFileName, workDirectoryName+"/")
    paramFileName, remain = readCommandFile(commandFileName)
    shutil.copy(paramFileName, workDirectoryName+"/")

    os.chdir(workDirectoryName)
    paramFileNameList = createMultipleParamFiles(paramFileName, sep, upLimit)
    pbsFileNameList = createScriptPBS(remain, paramFileNameList, commandFileName)
    submitJobs(pbsFileNameList, machine='server')


if __name__=='__main__':
    main()
