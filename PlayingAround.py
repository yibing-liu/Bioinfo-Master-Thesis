
import os
import re
import math
import pandas #read txt file as dataframe

#file path
rawCountsFilePath = os.getcwd() + "/Data/raw_counts"
segmentFilePath = os.getcwd() + "/Data/segmented"
perBinFilePath = os.getcwd() + "/Data/logr_cn_per_bin"
embryoFilePath = os.getcwd() + "/Output"

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

#get bin list
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
            embryoData = separateData[i] #initialise embryoData to proceed merge later
        embryoData = embryoData.merge (separateData[i])
    outputFileName = "TVEMB" + str(embryoNumber) + ".txt"
    embryoData.to_csv(outputFileName, sep="\t", mode="a") ###NEED TO SET OUTPUT-PATH
    return embryoData

##############################################################################################
#get EmbryodatafileName
def getEmbryoFile(embryoFileName):
    fileName = embryoFilePath + "/" + embryoFileName
    embryoFileDataFrame = pandas.read_csv(fileName, sep = "\t")
    return embryoFileDataFrame

#generate data frame for plotting density vs logR/bin (all together)
def getCombinedLogRData():
    combinedLogRData = pandas.DataFrame()
    separateData = dict()
    embryoNumber = 1
    cellNumber = 1
    combinedLogRData["Chr"] = getColumn(getRawFile(getRawCountsFileName(embryoNumber, cellNumber)), 0) #chromosome
    combinedLogRData["Start"] = getColumn(getRawFile(getRawCountsFileName(embryoNumber, cellNumber)), 1) #start
    combinedLogRData["End"] = getColumn(getRawFile(getRawCountsFileName(embryoNumber, cellNumber)), 2) #end
    combinedLogRData["Length"] = getColumn(getRawFile(getRawCountsFileName(embryoNumber, cellNumber)), 4) #length
    for i in range(1,48):
        for j in range(0,len(getEmbryoCellNumber(i))):
            combinedLogRData["TVEMB"+str(i)+"Cell"+str(getEmbryoCellNumber(i)[j])+"_logR"] = getColumn(getEmbryoFile("TVEMB" + str(i) + ".txt"), (4+3*j))
    combinedLogRData.to_csv("combinedLogRData.txt", sep="\t", mode="a")
    return combinedLogRData

#without #reads/total reads
def getSeparateCopyNumberData():
    copyNumber0Data = pandas.DataFrame(columns = ["chr", "start", "end", "CN", "CN_segment", "LogR", "LogR_segment"])
    copyNumber1Data = pandas.DataFrame(columns = ["chr", "start", "end", "CN", "CN_segment", "LogR", "LogR_segment"])
    copyNumber2Data = pandas.DataFrame(columns = ["chr", "start", "end", "CN", "CN_segment", "LogR", "LogR_segment"])
    copyNumber3Data = pandas.DataFrame(columns = ["chr", "start", "end", "CN", "CN_segment", "LogR", "LogR_segment"])
    for i in range(1,48): #first 5
        for j in range(0,len(getEmbryoCellNumber(i))):
            dataframe = getPerBinFile(getPerBinFileName(i, getEmbryoCellNumber(i)[j]))
            for n in range(0,len(dataframe)):
                currentRow = [dataframe.iloc[n]["chr"], dataframe.iloc[n]["start"], dataframe.iloc[n]["end"], dataframe.iloc[n]["CN"], dataframe.iloc[n]["CN_segment"], dataframe.iloc[n]["LogR"], dataframe.iloc[n]["LogR_segment"]]
                if (dataframe.iloc[n]["CN_segment"] == 0):
                    copyNumber0Data.loc[len(copyNumber0Data)] = currentRow
                if (dataframe.iloc[n]["CN_segment"] == 1):
                    copyNumber1Data.loc[len(copyNumber1Data)] = currentRow
                if (dataframe.iloc[n]["CN_segment"] == 2):
                    copyNumber2Data.loc[len(copyNumber2Data)] = currentRow
                if (dataframe.iloc[n]["CN_segment"] == 3):
                    copyNumber3Data.loc[len(copyNumber3Data)] = currentRow
    copyNumber0Data.to_csv("combinedCN0.txt", sep="\t", mode="a")
    copyNumber1Data.to_csv("combinedCN1.txt", sep="\t", mode="a")
    copyNumber2Data.to_csv("combinedCN2.txt", sep="\t", mode="a")
    copyNumber3Data.to_csv("combinedCN3.txt", sep="\t", mode="a")

# #reads/total reads
def getCurrentBinReads(embryoNumber, cellNumber, rowNumber):
    currentBinReads = getRawFile(getRawCountsFileName(embryoNumber, cellNumber)).loc[rowNumber-1][5]
    return currentBinReads

def getTotalBinReads(embryoNumber, cellNumber):
    totalBinReads = getColumn(getRawFile(getRawCountsFileName(embryoNumber, cellNumber)),5).sum()
    return totalBinReads

def getCorrectedSeparateCopyNumberData():
    correctedCopyNumber0Data = pandas.DataFrame(columns = ["chr", "start", "end", "CN", "CN_segment", "LogR", "LogR_segment", "LogR_corrected"])
    correctedCopyNumber1Data = pandas.DataFrame(columns = ["chr", "start", "end", "CN", "CN_segment", "LogR", "LogR_segment", "LogR_corrected"])
    correctedCopyNumber2Data = pandas.DataFrame(columns = ["chr", "start", "end", "CN", "CN_segment", "LogR", "LogR_segment", "LogR_corrected"])
    correctedCopyNumber3Data = pandas.DataFrame(columns = ["chr", "start", "end", "CN", "CN_segment", "LogR", "LogR_segment", "LogR_corrected"])
    for i in range(1,48):
        for j in range(0,len(getEmbryoCellNumber(i))):
            dataframe = getPerBinFile(getPerBinFileName(i, getEmbryoCellNumber(i)[j]))
            for n in range(0,len(dataframe)):
                LogR_corrected = dataframe.iloc[n]["LogR"]*getCurrentBinReads(i, getEmbryoCellNumber(i)[j], (n+1))/getTotalBinReads(i, getEmbryoCellNumber(i)[j])
                currentRow = [dataframe.iloc[n]["chr"], dataframe.iloc[n]["start"], dataframe.iloc[n]["end"], dataframe.iloc[n]["CN"], dataframe.iloc[n]["CN_segment"], dataframe.iloc[n]["LogR"], dataframe.iloc[n]["LogR_segment"], LogR_corrected]
                if (dataframe.iloc[n]["CN_segment"] == 0):
                    correctedCopyNumber0Data.loc[len(correctedCopyNumber0Data)] = currentRow
                if (dataframe.iloc[n]["CN_segment"] == 1):
                    correctedCopyNumber1Data.loc[len(correctedCopyNumber1Data)] = currentRow
                if (dataframe.iloc[n]["CN_segment"] == 2):
                    correctedCopyNumber2Data.loc[len(correctedCopyNumber2Data)] = currentRow
                if (dataframe.iloc[n]["CN_segment"] == 3):
                    correctedCopyNumber3Data.loc[len(correctedCopyNumber3Data)] = currentRow
    correctedCopyNumber0Data.to_csv("combinedCN0_corrected.txt", sep="\t", mode="a")
    correctedCopyNumber1Data.to_csv("combinedCN1_corrected.txt", sep="\t", mode="a")
    correctedCopyNumber2Data.to_csv("combinedCN2_corrected.txt", sep="\t", mode="a")
    correctedCopyNumber3Data.to_csv("combinedCN3_corrected.txt", sep="\t", mode="a")

#################################################################################################
#get probability from a poisson distribution, given lambda
def getPoissonProbability(l, k):
    probability = (l**k) * math.exp(-l)/math.gamma(k+1)
    return probability

#generate file logR -->probability(cn0, cn1, cn2, cn3)
###############################################################################################
def getProbabilityLogRData():
    for i in range(1,48):
        for j in range(0,len(getEmbryoCellNumber(i))):
            probabilityLogRData = pandas.DataFrame(columns = ["chr", "start", "end", "LogR_corrected", "Pcn0", "Pcn1", "Pcn2", "Pcn3"])
            dataframe = getPerBinFile(getPerBinFileName(i, getEmbryoCellNumber(i)[j]))
            for n in range(0,len(dataframe)):
                LogR_corrected = dataframe.iloc[n]["LogR"]*getCurrentBinReads(i, getEmbryoCellNumber(i)[j], (n+1))/getTotalBinReads(i, getEmbryoCellNumber(i)[j])
                Pcn0 = getPoissonProbability(-0.0002759783, LogR_corrected) ##correct lambda value needed
                Pcn1 = getPoissonProbability(0.9999388, LogR_corrected+1)   ##correct lambda value needed
                Pcn2 = getPoissonProbability(1.999995, LogR_corrected+2)    ##correct lambda value needed
                Pcn3 = getPoissonProbability(0.7042147, LogR_corrected+3)   ##correct lambda value needed
                currentRow = [dataframe.iloc[n]["chr"], dataframe.iloc[n]["start"], dataframe.iloc[n]["end"], LogR_corrected, Pcn0, Pcn1, Pcn2, Pcn3]
                print(currentRow)
                probabilityLogRData.loc[len(probabilityLogRData)] = currentRow
                outputFileName = "TVEMB" + str(i) + "_" + getEmbryoCellNumber(i)[j] + "_probability.txt"
            probabilityLogRData.to_csv(outputFileName, sep="\t", mode="a")


###############################################################################################
###############################################################################################
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

#############################################################################################
#CONTROL
def getControlNumber(fileType):
    cellNumber = []
    if fileType == "control":
        for file in os.listdir(rawCountsFilePath):
            fileName = os.path.basename(file)
            if fileName.startswith("CONTROL"):
                currentCellNumber = (fileName.split("CONTROL_")[1]).split('.')[0]
                cellNumber.append(currentCellNumber)
                cellNumber.sort(key = int)
        return cellNumber
    if fileType == "empty":
        cellNumber = []
        return cellNumber
    if fileType == "NC":
        for file in os.listdir(rawCountsFilePath):
            fileName = os.path.basename(file)
            if fileName.startswith("NC"):
                currentCellNumber = (fileName.split("NC_")[1]).split('.')[0]
                cellNumber.append(currentCellNumber)
                cellNumber.sort(key = int)
        return cellNumber
    if fileType == "PCMC":
        for file in os.listdir(rawCountsFilePath):
            fileName = os.path.basename(file)
            if fileName.startswith("PC_MC"):
                currentCellNumber = (fileName.split("PC_MC_")[1]).split('.')[0]
                cellNumber.append(currentCellNumber)
                cellNumber.sort(key = int)
        return cellNumber
    if fileType == "PCSC":
        for file in os.listdir(rawCountsFilePath):
            fileName = os.path.basename(file)
            if fileName.startswith("PC_SC"):
                currentCellNumber = (fileName.split("PC_SC_")[1]).split('.')[0]
                cellNumber.append(currentCellNumber)
                cellNumber.sort(key = int)
        return cellNumber

def pasteColumn(columnNumber, rawFileName, perBinFileName):
    separateData = dict()
    separateData[columnNumber] = pandas.DataFrame()
    separateData[columnNumber]["Chr"] = getColumn(getRawFile(rawFileName), 0) #chromosome
    separateData[columnNumber]["Start"] = getColumn(getRawFile(rawFileName), 1) #start
    separateData[columnNumber]["End"] = getColumn(getRawFile(rawFileName), 2) #end
    separateData[columnNumber]["Length"] = getColumn(getRawFile(rawFileName), 4) #length
    separateData[columnNumber]["cell"+str(columnNumber)+"_binReads"] = getColumn(getRawFile(rawFileName), 5) #binReads
    separateData[columnNumber]["cell"+str(columnNumber)+"_CN"] = getColumn(getPerBinFile(perBinFileName), 3) #cn
    separateData[columnNumber]["cell"+str(columnNumber)+"_logR"] = getColumn(getPerBinFile(perBinFileName), 5) #logR
    return separateData[columnNumber]

def getControlData(fileType):
    controlData = pandas.DataFrame()
    separateData = dict()
    if fileType == "control":
        for i in getControlNumber(fileType):
            rawFileName = "CONTROL_" + str(i) + ".250K.101.sorted.count-t"
            perBinFileName = "CONTROL_" + str(i) +".250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.txt"
            if i == getControlNumber(fileType)[0]:
                controlData = pasteColumn(i, rawFileName, perBinFileName)
            controlData = controlData.merge (pasteColumn(i, rawFileName, perBinFileName))
        outputFileName = "control" + ".txt"
        controlData.to_csv(outputFileName, sep='\t', mode='a') ###NEED TO SET OUTPUT-PATH
        return controlData
    if fileType == "empty":
        rawFileName = "empty_.250K.101.sorted.count-t"
        perBinFileName = "empty_.250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.txt"
        controlData = pasteColumn(0, rawFileName, perBinFileName)
        outputFileName = "empty" + ".txt"
        controlData.to_csv(outputFileName, sep='\t', mode='a') ###NEED TO SET OUTPUT-PATH
        return controlData
    if fileType == "NC":
        for i in getControlNumber(fileType):
            rawFileName = "NC_" + str(i) + ".250K.101.sorted.count-t"
            perBinFileName = "NC_" + str(i) +".250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.txt"
            if i == getControlNumber(fileType)[0]:
                controlData = pasteColumn(i, rawFileName, perBinFileName)
            controlData = controlData.merge (pasteColumn(i, rawFileName, perBinFileName))
        outputFileName = "NC" + ".txt"
        controlData.to_csv(outputFileName, sep='\t', mode='a') ###NEED TO SET OUTPUT-PATH
        return controlData
    if fileType == "PCMC":
        for i in getControlNumber(fileType):
            rawFileName = "PC_MC_" + str(i) + ".250K.101.sorted.count-t"
            perBinFileName = "PC_MC_" + str(i) +".250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.txt"
            if i == getControlNumber(fileType)[0]:
                controlData = pasteColumn(i, rawFileName, perBinFileName)
            controlData = controlData.merge (pasteColumn(i, rawFileName, perBinFileName))
        outputFileName = "PC_MC" + ".txt"
        controlData.to_csv(outputFileName, sep='\t', mode='a') ###NEED TO SET OUTPUT-PATH
        return controlData
    if fileType == "PCSC":
        for i in getControlNumber(fileType):
            rawFileName = "PC_SC_" + str(i) + ".250K.101.sorted.count-t"
            perBinFileName = "PC_SC_" + str(i) +".250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.txt"
            if i == getControlNumber(fileType)[0]:
                controlData = pasteColumn(i, rawFileName, perBinFileName)
            controlData = controlData.merge (pasteColumn(i, rawFileName, perBinFileName))
        outputFileName = "PC_SC" + ".txt"
        controlData.to_csv(outputFileName, sep='\t', mode='a') ###NEED TO SET OUTPUT-PATH
        return controlData


#############################################################################################
#OPERATOR
#for i in range(1,48):
#    getEmbryoData(i)

#getCombinedLogRData()

#getSeparateCopyNumberData()
#getCorrectedSeparateCopyNumberData()
getProbabilityLogRData()


#fileType = ["control", "empty", "NC", "PCMC", "PCSC"]
#for i in fileType:
#    getControlData(i)
