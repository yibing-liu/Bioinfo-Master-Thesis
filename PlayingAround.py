
import os
import pandas #read txt file as dataframe

#file path
rawCountsFilePath = os.getcwd()+"/Data/raw_counts"
segmentFilePath = os.getcwd()+"/Data/segmented"
perBinFilePath = os.getcwd()+"/Data/logr_cn_per_bin"

#generate raw count file fileName
def getRawCountsFileName(embryoNumber, i):
    rawCountsFileName = rawCountsFilePath + "TVEMB" + embryoNumber + "_" + i +".250K.101.sorted.count-t"
    return rawCountsFileName

#generate segmented file fileName
def getSegmentFileName(embryoNumber, i):
    segmentFileName = segmentFilePath + "TVEMB" + embryoNumber + "_" + i +".250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.breakpoints.txt"
    return segmentFileName

#generate per bin file fileName
def getPerBinFileName(embryoNumber, i):
    perBinFileName = perBinFilePath + "TVEMB" + embryoNumber + "_" + i +".250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.txt"
    return perBinFileName

#read in raw count file as dataframe
def getRawFile(rawFileName):
    rawCountsDataFrame = pandas.read_csv(rawFileName, sep = "\t", header = None)
    rawCountsDataFrame.columns = ["chr", "start", "end", "binSize", "length", "binReads", "systemControl"]
    return rawCountsDataFrame

#read in segment file as dataframe
def getSegmentFile(segmentFileName):
    segmentDataFrame = pandas.read_csv(segmentFileName, sep = "\t")
    return segmentDataFrame

#read in logR cn /bin
def getPerBinFile(perBinFileName):
    perBinDataFrame = pandas.read_csv(perBinFileName, sep = "\t")
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

#getColumn(fileName, 1)
def getColumn(fileName, columnNumber):
    columnName = list(file)[columnNumber]
    column = file[columnName]
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
    n = getEmbryoCellNumber(embryoNumber)[0]
    columnBinChr = getColumn(getRawFile(getRawCountsFileName(embryoNumber, n)), 1) #chromosome
    columnBinStart = getColumn(getRawFile(getRawCountsFileName(embryoNumber, n)), 2) #start
    columnBinEnd = getColumn(getRawFile(getRawCountsFileName(embryoNumber, n)), 3) #end
    columnBinLength = getColumn(getRawFile(getRawCountsFileName(embryoNumber, n)), 5) #length
    embryoData.append(columnBinChr, columnBinStart, columnBinEnd, columnBinLength)
    for i in getEmbryoCellNumber(embryoNumber):
        columnBinReads = getColumn(getRawFile(getRawCountsFileName(embryoNumber, i)), 6) #binReads
        columnBinReads.columns = ["cell"+i+"_binReads"]
        embryoData.append(columnBinReads)
        columnCN = getColumn(getPerBinFile(getPerBinFileName(embryoNumber, i)), 4) #cn
        columnCN.columns = ["cell"+i+"_CN"]
        embryoData.append(columnCN)
        columnLogR = getColumn(getPerBinFile(getPerBinFileName(embryoNumber, i)), 6) #logR
        columnLogR.columns = ["cell"+i+"_logR"]
        embryoData.append(columnLogR)
    #outputFileName = "embryo"+embryoNumber
    #embryoData.to_csv(outputFileName, sep='\t', mode='a')
    print(embryoData)
#input
def input():
    embryoNumber = "4"#input('Embryo Number:')  ### BUG:
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
tester = getEmbryoCellNumber(1)







print(tester)
