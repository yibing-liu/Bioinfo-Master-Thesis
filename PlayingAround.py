
import os

#files in raw_counts folder
rawCountsFile = []
path = os.getcwd()+"/data/raw_counts"
print(path)
for i in os.listdir(path):
    fileName = os.path.basename(i)
    print(fileName)
    #files.append(open())
    #os.path.splitext(path)
