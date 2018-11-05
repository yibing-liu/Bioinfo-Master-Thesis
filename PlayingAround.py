
import os
import re
import pandas #read txt file as dataframe

#file path
rawCountsFilePath = os.getcwd()+"/Data/raw_counts"
segmentFilePath = os.getcwd()+"/Data/segmented"
perBinFilePath = os.getcwd()+"/Data/logr_cn_per_bin"

#generate raw count file fileName
def getRawCountsFileName(embryoNumber, i):
    rawCountsFileName = "TVEMB" + str(embryoNumber) + "_" + str(i) +".250K.101.sorted.count-t"
    return rawCountsFileName

#generate segmented file fileName
def getSegmentFileName(embryoNumber, i):
    segmentFileName = "TVEMB" + str(embryoNumber) + "_" + str(i) +".250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.breakpoints.txt"
    return segmentFileName

#generate per bin file fileName
def getPerBinFileName(embryoNumber, i):
    perBinFileName = "TVEMB" + str(embryoNumber) + "_" + str(i) +".250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.txt"
    return perBinFileName

#read in raw count file as dataframe
def getRawFile(rawFileName):
    fileName = rawCountsFilePath + "/" + rawFileName
    rawCountsDataFrame = pandas.read_csv(fileName, sep = "\t", header = None)
    rawCountsDataFrame.columns = ["chr", "start", "end", "binSize", "length", "binReads", "systemControl"]
    return rawCountsDataFrame

#read in segment file as dataframe
def getSegmentFile(segmentFileName):
    fileName = segmentFilePath + "/" + segmentFileName
    segmentDataFrame = pandas.read_csv(fileName, sep = "\t")
    return segmentDataFrame

#read in logR cn /bin
def getPerBinFile(perBinFileName):
    fileName = perBinFilePath + "/" + perBinFileName
    perBinDataFrame = pandas.read_csv(fileName, sep = "\t")
    return perBinDataFrame

#get files of given embryo
def getFile(inputFileName):
    if bool("sorted.count" in inputFileName):
        path = rawCountsFilePath
        for file in os.listdir(path):
            fileName = os.path.basename(file)
            if inputFileName == fileName:
                outputFile = getRawFile(inputFileName)
    if bool("gc_corrected.segments.copynumber.breakpoints" in inputFileName):
        path = segmentFilePath
        for file in os.listdir(path):
            fileName = os.path.basename(file)
            if inputFileName == fileName:
                outputFile = getSegmentFile(inputFileName)
    if bool("sorted.count.gamma_15.gc_corrected.segments.copynumber" in inputFileName):
        path = perBinFilePath
        for file in os.listdir(path):
            fileName = os.path.basename(file)
            if inputFileName == fileName:
                outputFile = getPerBinFile(inputFileName)
    return outputFile

#get bin list  ##### BUG:
def getBinList(embryoNumber):
    if (getRawFile(embryoNumber, 1) != None):
        fileName = getRawFile(embryoNumber, 1) #certain embryo's cell #1 doesnt exist
    else:
        fileName = getRawFile(embryoNumber, 2)
    binList = getRawFile(rawFileName)[1, 2, 3]
    return binList

#getColumn(fileName, 1)
def getColumn(file, columnNumber):
    column = file[list(file)[columnNumber]]
    return column

#how many cells are linked to selected embryo
def getEmbryoCellNumber(embryoNumber):
    embryoFileName = "TVEMB" + str(embryoNumber) + "_"
    cellNumber = []
    for file in os.listdir(rawCountsFilePath):
        fileName = os.path.basename(file)
        if fileName.startswith(embryoFileName):
            currentCellNumber = (fileName.split(embryoFileName)[1]).split('.')[0]
            cellNumber.append(currentCellNumber)
            cellNumber.sort(key = int)
    return cellNumber


#output file per embryo
def getEmbryoData(embryoNumber):
    embryoData = pandas.DataFrame()
    separateData = dict()
    for i in getEmbryoCellNumber(embryoNumber):
        separateData[i] = pandas.DataFrame()
        separateData[i]["Chr"] = getColumn(getRawFile(getRawCountsFileName(embryoNumber, i)), 0) #chromosome
        separateData[i]["Start"] = getColumn(getRawFile(getRawCountsFileName(embryoNumber, i)), 1) #start
        separateData[i]["End"] = getColumn(getRawFile(getRawCountsFileName(embryoNumber, i)), 2) #end
        separateData[i]["Length"] = getColumn(getRawFile(getRawCountsFileName(embryoNumber, i)), 4) #length
        separateData[i]["cell"+i+"_binReads"] = getColumn(getRawFile(getRawCountsFileName(embryoNumber, i)), 5) #binReads
        separateData[i]["cell"+i+"_CN"] = getColumn(getPerBinFile(getPerBinFileName(embryoNumber, i)), 3) #cn
        separateData[i]["cell"+i+"_logR"] = getColumn(getPerBinFile(getPerBinFileName(embryoNumber, i)), 5) #logR
        if i == getEmbryoCellNumber(embryoNumber)[0]:
            embryoData = separateData[i]
    embryoData = embryoData.merge (separateData[i])
    outputFileName = "output_embryo" + str(embryoNumber) + ".txt"
    embryoData.to_csv(outputFileName, sep='\t', mode='a') ###NEED TO SET OUTPUT-PATH
    return embryoData

#input fileName
def inputparameter():
    embryoNumber = input('Embryo Number:')
    cellNumber = input('cell Number:')
    fileType = input('File Type:')
    if fileType == "raw_counts":
        inputfileName = getRawCountsFileName(embryoNumber, cellNumber)
    if fileType == "segmented":
        inputfileName = getSegmentFileName(embryoNumber, cellNumber)
    if fileType == "per bin":
        inputfileName = getPerBinFileName(embryoNumber, cellNumber)
    return inputfileName

##TESTER
print(getEmbryoData(9))
