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
#outputFile_cn2_bin = open('cn2_bin.txt', 'a+')
for i in range(0,len(mapD_selected)):
    if 'TVEMB' in mapD_selected.iloc[i]['sample']:
        embryoNumber = mapD_selected.iloc[i]['sample'].split("_")[0].split("TVEMB")[1]
        cellNumber = mapD_selected.iloc[i]['sample'].split("_")[1]
        for line in enumerate(open(perBinFilePath + getPerBinFileName(embryoNumber, cellNumber))):
            outputFileName_cn1 = 'cn1.txt'
            outputFileName_cn2 = 'cn2.txt'
            outputFileName_cn3 = 'cn3.txt'
            outputFile_cn1_bin = open(outputFileName_cn1, 'a+')
            outputFile_cn2_bin = open(outputFileName_cn2, 'a+')
            outputFile_cn3_bin = open(outputFileName_cn3, 'a+')
            if line[1].split()[4] == '1':
                outputFile_cn1_bin.write(line[1].split()[0]+'\t'+line[1].split()[1]+'\t'+line[1].split()[2]+'\t'+line[1].split()[5]+'\n')
            elif line[1].split()[4] == '2':
                outputFile_cn2_bin.write(line[1].split()[0]+'\t'+line[1].split()[1]+'\t'+line[1].split()[2]+'\t'+line[1].split()[5]+'\n')
            elif line[1].split()[4] == '3':
                outputFile_cn3_bin.write(line[1].split()[0]+'\t'+line[1].split()[1]+'\t'+line[1].split()[2]+'\t'+line[1].split()[5]+'\n')



#print(open("tester.txt"))
