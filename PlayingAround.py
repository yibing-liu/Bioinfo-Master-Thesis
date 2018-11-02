
import os
import pandas #read txt file as dataframe

#file path
rawCountsFilePath = os.getcwd()+"/Data/raw_counts"
segmentFilePath = os.getcwd()+"/Data/segmented"
perBinFilePath = os.getcwd()+"/Data/logr_cn_per_bin"

#generate file fileName
def getRawCountsFileName(embryoNumber, i):
    rawCountsFileName = "TVEMB" + embryoNumber + "_" + i +".250K.101.sorted.count-t"
    return rawCountsFileName


#get raw_counts files of given embryo
def getFile(FileName):
    if fileName.find("sorted.count"):
        path = rawCountsFilePath
    if fileName.find("gc_corrected.segments.copynumber.breakpoints"):
        path = segmentFilePath
    if fileName.find("sorted.count.gamma_15.gc_corrected.segments.copynumber"):
        path = perBinFilePath
    for file in os.listdir(path):
        fileName = os.path.basename(file)
        if fileName == rawCountsFileName:
            return file

#read in raw count file as dataframe
def getRawFile(rawFileName):
    rawCountsDataFrame = pandas.read_csv(rawFileName, header = None)
    dataframe.columns = ["chr", "start", "end", "binSize", "length", "binReads", "systemControl"]
    return rawCountsDataFrame

#read in segment file as dataframe
def getSegmentFile(segmentFileName):
    segmentDataFrame = pandas.read_csv(segmentFileName, header = True)
    return segmentDataFrame

#read in logR cn /bin
def getPerBinFile(perBinFileName):
    perBinDataFrame = pandas.read_csv(perBinFileName, header = True)
    return perBinFileName

#paste interested columns in file
def















#files in raw_counts folder
rawCountsFile = []
path = os.getcwd()+"/Data/raw_counts"
print(path)
for i in os.listdir(path):
    fileName = os.path.basename(i)
    print(fileName)
    files.append(open())
    os.path.splitext(path)
