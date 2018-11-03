
import os
import pandas #read txt file as dataframe

#file path
rawCountsFilePath = os.getcwd()+"/Data/raw_counts"
segmentFilePath = os.getcwd()+"/Data/segmented"
perBinFilePath = os.getcwd()+"/Data/logr_cn_per_bin"

#generate raw count file fileName
def getRawCountsFileName(embryoNumber, i):
    rawCountsFileName = "TVEMB" + embryoNumber + "_" + i +".250K.101.sorted.count-t"
    return rawCountsFileName

#generate segmented file fileName
def getSegmentFileName(embryoNumber, i):
    segmentFileName = "TVEMB" + embryoNumber + "_" + i +".250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.breakpoints.txt"
    return segmentFileName

#generate per bin file fileName
def getPerBinFileName(embryoNumber, i):
    perBinFileName = "TVEMB" + embryoNumber + "_" + i +".250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.txt"
    return perBinFileName

#read in raw count file as dataframe
def getRawFile(rawFileName):
    rawCountsDataFrame = pandas.read_csv(rawFileName, header = None)
    dataframe.columns = ["chr", "start", "end", "binSize", "length", "binReads", "systemControl"]
    return rawCountsDataFrame

#read in segment file as dataframe
def getSegmentFile(segmentFileName):
    segmentDataFrame = pandas.read_csv(segmentFileName, header = None) #True) BUG
    return segmentDataFrame

#read in logR cn /bin
def getPerBinFile(perBinFileName):
    perBinDataFrame = pandas.read_csv(perBinFileName, header = None) #True) BUG
    return perBinFileName

#get raw_counts files of given embryo
def getFile(inputFileName):
    if inputFileName.find("sorted.count"):
        path = rawCountsFilePath + "/" + inputFileName
        for file in os.listdir(path):
            fileName = os.path.basename(file)
            if inputFileName == fileName:
                getRawFile(inputFileName)
    if fileName.find("gc_corrected.segments.copynumber.breakpoints"):
        path = segmentFilePath + "/" + inputFileName
        for file in os.listdir(path):
            fileName = os.path.basename(file)
            if inputFileName == fileName:
                getSegmentFile(inputFileName)
    if fileName.find("sorted.count.gamma_15.gc_corrected.segments.copynumber"):
        path = perBinFilePath + "/" + inputFileName
        for file in os.listdir(path):
            fileName = os.path.basename(file)
            if inputFileName == fileName:
                getPerBinFile(inputFileName)

#get bin list
def getBinList(embryoNumber):
    if (getRawFile(embryoNumber, 1) != None):
        fileName = getRawFile(embryoNumber, 1) #certain embryo's cell #1 doesnt exist
    else:
        fileName = getRawFile(embryoNumber, 2)
    binList = getRawFile(rawFileName)[1, 2, 3]
    return binList

#paste interested columns in file
def pasteColumn(fileName, columnNumber): #, outputFileName
    file = getFile(fileName)
    outputFile = pandas.dataframe(file, columns = [columnNumber])
    print(outputFile)

#inputFileName
def input():
    embryoNumber = "1"#input('Embryo Number:')  ### BUG:
    cellNumber = "1"#input('Cell Number:')  ##BUG
    fileType = "raw_counts"#input('File Type:') ##BUG
    if fileType == "raw_counts":
        fileName = getRawCountsFileName(embryoNumber, cellNumber)
    if fileType == "segmented":
        fileName = getSegmentFileName(embryoNumber, cellNumber)
    if fileType == "per bin":
        fileName = getPerBinFileName(embryoNumber, cellNumber)
    return fileName

##TESTER
fileName = input()

print(fileName)
print("TVEMB1_1.250K.101.sorted.count-t")
print(rawCountsFilePath)
rawCountsFilePath = rawCountsFilePath + "/" + fileName
print(rawCountsFilePath)
file = pandas.read_csv(rawCountsFilePath) ##NEED TO FIX SEP
print(file)
#pasteColumn(fileName, 1)
