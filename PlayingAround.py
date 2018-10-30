
import os

#generate file fileName
def getRawCountsFileName(a, i):
    rawCountsFileName = "TVEMB" + a + "_" + i +".250K.101.sorted.count-t"
    return rawCountsFileName

#get raw_counts files of given embryo
def getRawCountsFile(rawCountsFileName):
    for file in os.listdir(path):
        fileName = os.path.basename(file)
        if fileName == rawCountsFileName:
            return file

#select interested columns in file
def getColumn(rawCountsFile):
    













#files in raw_counts folder
rawCountsFile = []
path = os.getcwd()+"/Data/raw_counts"
print(path)
for i in os.listdir(path):
    fileName = os.path.basename(i)
    print(fileName)
    files.append(open())
    os.path.splitext(path)
