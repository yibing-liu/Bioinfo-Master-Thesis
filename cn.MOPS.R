source("https://bioconductor.org/biocLite.R")
biocLite("cn.mops")
library(cn.mops)
getwd()
tvemb1 = read.table("/Users/Yibing/Bioinfo-Master-Thesis/Output/TVEMB1.txt", sep = "\t", header = TRUE)
#tvemb2 = read.table("/Users/Yibing/Bioinfo-Master-Thesis/Output/TVEMB2.txt", sep = "\t", header = TRUE)
rowName = read.table("/Users/Yibing/Bioinfo-Master-Thesis/Output/rowName.txt", header = TRUE)
tvemb1 = cbind(tvemb1, rowName)
chrName = list("1","2","3","4","5","6","7","8","9","10",
            "11","12","13","14","15","16","17","18","19","20","21","22","X")
for(i in chrName){
  subsetName = paste("tvemb1.chr", i, sep = "")
  data = subset(tvemb1, tvemb1$Chr==i)
  assign(subsetName, data)
  }


tvemb1.binreads.chr1 = data.matrix(tvemb1.chr1[c(6, 9)])
rownames(tvemb1.binreads.chr1) = tvemb1.chr1$Chr_Start_End_Length
resCNMOPSX = cn.mops(tvemb1.binreads.chr1)
segplot(resCNMOPSX)

#embro46
tvemb46 = read.table("/Users/Yibing/Bioinfo-Master-Thesis/Output/TVEMB46.txt", sep = "\t", header = TRUE)
tvemb46 = cbind(tvemb46, rowName)
for(i in chrName){
  subsetName = paste("tvemb46.chr", i, sep = "")
  data = subset(tvemb46, tvemb46$Chr==i)
  assign(subsetName, data)
}
tvemb46.binreads.chr1 = data.matrix(tvemb46.chr1[c(seq(6, 77, by=3))])
rownames(tvemb46.binreads.chr1) = tvemb46.chr1$Chr_Start_End_Length
resCNMOPSX = cn.mops(tvemb46.binreads.chr1)
segplot(resCNMOPSX)


#embro40
tvemb40 = read.table("/Users/Yibing/Bioinfo-Master-Thesis/Output/TVEMB40.txt", sep = "\t", header = TRUE)
tvemb40 = cbind(tvemb40, rowName)
for(i in chrName){
  subsetName = paste("tvemb40.chr", i, sep = "")
  data = subset(tvemb40, tvemb40$Chr==i)
  assign(subsetName, data)
}
tvemb40.binreads.chr3 = data.matrix(tvemb40.chr3[c(seq(24, 48, by=3))])
rownames(tvemb40.binreads.chr3) = tvemb40.chr3$Chr_Start_End_Length
resCNMOPSX = cn.mops(tvemb40.binreads.chr3)
segplot(resCNMOPSX)


testx = cn.mops(resCNMOPSX)
resCNMOPS <- calcIntegerCopyNumbers(resCNMOPSX)
plot(resCNMOPS, which = 1)


plot(resCNMOPSX)

resCNMOPSX = cn.mops(X)

plot(resCNMOPS, which = 1)

tvemb46.binreads.chr19 = data.matrix(tvemb46.chr19[c(9, 12, 27, 33, 54)])
rownames(tvemb46.binreads.chr19) = tvemb46.chr19$Chr_Start_End_Length
resCNMOPSX = cn.mops(tvemb46.binreads.chr19)
segplot(resCNMOPSX)

tvemb46.binreads.chr3 = data.matrix(tvemb46.chr3[c(9, 12, 27, 33, 54)])
rownames(tvemb46.binreads.chr3) = tvemb46.chr3$Chr_Start_End_Length
resCNMOPSX = cn.mops(tvemb46.binreads.chr3)
segplot(resCNMOPSX)



#embro47
tvemb47 = read.table("/Users/Yibing/Bioinfo-Master-Thesis/Output/TVEMB47.txt", sep = "\t", header = TRUE)
tvemb47 = cbind(tvemb47, rowName)
for(i in chrName){
  subsetName = paste("tvemb47.chr", i, sep = "")
  data = subset(tvemb47, tvemb47$Chr==i)
  assign(subsetName, data)
}
tvemb47.binreads.chr1 = data.matrix(tvemb47.chr1[c(seq(51, 95, by=3))])
rownames(tvemb47.binreads.chr1) = tvemb47.chr1$Chr_Start_End_Length
resCNMOPSX = cn.mops(tvemb47.binreads.chr1)
segplot(resCNMOPSX)

#combine embro 1+46+47
combine.binreads.chr1 = cbind(tvemb1.chr1[c(6)], tvemb46.chr1[c(seq(6, 77, by=3))], tvemb47.chr1[c(seq(6, 95, by=3))])
combine.binreads.chr1 = data.matrix(combine.binreads.chr1)
rownames(combine.binreads.chr1) = tvemb47.chr1$Chr_Start_End_Length
colnames(combine.binreads.chr1) = c(1:55)
resCNMOPSX = cn.mops(combine.binreads.chr1)
segplot(resCNMOPSX)

combine.binreads.chr9 = cbind(tvemb1.chr9[c(6)], tvemb46.chr9[c(seq(6, 77, by=3))], tvemb47.chr9[c(seq(6, 95, by=3))])
combine.binreads.chr9 = data.matrix(combine.binreads.chr9)
rownames(combine.binreads.chr9) = tvemb47.chr9$Chr_Start_End_Length
colnames(combine.binreads.chr9) = c(1:56)
resCNMOPSX = cn.mops(combine.binreads.chr9)
segplot(resCNMOPSX)

testx = cn.mops(X)
resCNMOPS <- calcIntegerCopyNumbers(testx)
plot(resCNMOPS, which = 1)

