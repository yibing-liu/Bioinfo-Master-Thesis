import os
import re
import numpy
import pandas

#file path
perBinFilePath = os.getcwd() + '/Data/logr_cn_per_bin/'

#generate per bin file fileName
def getPerBinFileName(embryoNumber, i):
    perBinFileName = 'TVEMB' + str(embryoNumber) + '_' + str(i) + '.250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.txt'
    return perBinFileName

def getPerBinFile(perBinFileName):
    fileName = perBinFilePath + "/" + perBinFileName
    perBinDataFrame = pandas.read_csv(fileName, sep = "\t")
    return perBinDataFrame

def getColumn(file, columnNumber):
    column = file[list(file)[columnNumber]]
    return column

mapDFileName = os.getcwd()+'/mapD_selected.csv'
mapD_selected = pandas.read_csv(mapDFileName)


#get all bins with cn=2
#outputFile_cn2_bin = open('cn2_bin.txt', 'w+')
for i in range(0,len(mapD_selected)):
    sampleName = mapD_selected.iloc[i]['sample']
    split = sampleName.split("_")
    if 'TVEMB' in split[0]:
        embryoNumber = split[0].split("TVEMB")[1]
        cellNumber = split[1]
        filepath = perBinFilePath + getPerBinFileName(embryoNumber, cellNumber)
        median = numpy.median(getColumn(getPerBinFile(getPerBinFileName(embryoNumber, cellNumber)), 5))
        print(median)
        for line in enumerate(open(filepath)):
            if line[1].split()[4] == '2':
                outputFileName = 'cn2_TVEMB'+embryoNumber+'_'+cellNumber+'.txt'
                #outputFile_cn2_bin = open(outputFileName, 'w+')
                #outputFile_cn2_bin.write(line[1])
#outputFile.close()



for line in enumerate(open(perBinFilePath + getPerBinFileName(1, 1))):
    if line[1].split()[4] == '2':
        print(line)
