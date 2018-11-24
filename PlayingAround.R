install.packages("dplyr")
library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)
#Poisson
install.packages("fitdistrplus")
library(fitdistrplus)


getwd()

#plot all embryo cells logR+cn
for (i in 1:47){
  fileName = paste(paste("/Users/Yibing/Bioinfo-Master-Thesis/Output/TVEMB", i, sep = ""), "txt", sep = ".")
  data = read.csv(fileName, sep = "\t", header = TRUE)
  data$position = (data$End - data$Start)/2 + data$Start
  
  chromosomeNumber = c(seq(1,22), "X")
  columnNumber_logR = c(seq(8,ncol(data)-1,3), ncol(data))
  columnNumber_cn = c(seq(7,ncol(data)-2,3), ncol(data))
  plotDataCombined_logR = data.frame()
  plotDataCombined_cn = data.frame()
  
  for (j in chromosomeNumber){
    plotData_logR = data.frame(filter(data, (data$Chr==j)))[columnNumber_logR]
    plotData_cn = data.frame(filter(data, (data$Chr==j)))[columnNumber_cn]
    logR = melt(plotData_logR, "position")
    cn = melt(plotData_cn, "position")
    logR$Chr = j
    cn$Chr = j
    plotDataCombined_logR = rbind(plotDataCombined_logR, logR)
    plotDataCombined_cn = rbind(plotDataCombined_cn, cn)
  }
  
  plotDataName_logR = paste(paste("TVEMB", i, sep = ""), "_logR", sep = "")
  plotDataName_cn = paste(paste("TVEMB", i, sep = ""), "_cn", sep = "")
  
  pdf(plotDataName_logR, height = 16, width = 12)
  plot = ggplot(plotDataCombined_logR, aes(position, value, group=variable, color = variable)) +
    geom_point(size = 0.2) + 
    facet_wrap(~Chr, scales = "free_x", ncol = 3) +
    labs(x = "Position", y = "logR", color = NULL) +
    theme(legend.position = "bottom")
  print(plot)
  dev.off()
  
  png(plotDataName_cn, height = 16, width = 12)
  plot = ggplot(plotDataCombined_cn, aes(position, value, group=variable, color = variable)) +
    geom_point(size = 0.2) + 
    facet_wrap(~Chr, scales = "free_x", ncol = 3) +
    labs(x = "Position", y = "cn", color = NULL) +
    theme(legend.position = "bottom")
  print(plot)
  dev.off()
}


#### order chromosome!!!!
#### x y axis scale!!!!
####graph region enlarge!!!!
############################################################################
#density vs logR
#all cells together
dataCombined = read.csv("/Users/Yibing/Bioinfo-Master-Thesis/Output/combinedLogRData.txt", sep = "\t", header = TRUE)
dataCombined$position = (dataCombined$End - dataCombined$Start)/2 + dataCombined$Start
#dim(dataCombined)
dataCombined = dataCombined[c(6:322)]
logRCombined = melt(dataCombined, "position")
ggplot(logRCombined, aes(value)) + geom_density() + labs(x = "logR", y = "Density", color = NULL) 

poisson_logR = fitdist(logRCombined$value, 'pois', method = 'mme')   ?????methods


##########################################################################################
#density vs logR
#separate copy number(0, 1, 2, 3) vs logR (without #reads/total reads)
dataCopyNumber0 = read.csv("/Users/Yibing/Bioinfo-Master-Thesis/Output/combinedCN0.txt", sep = "\t", header = TRUE)
dataCopyNumber1 = read.csv("/Users/Yibing/Bioinfo-Master-Thesis/Output/combinedCN1.txt", sep = "\t", header = TRUE)
dataCopyNumber2 = read.csv("/Users/Yibing/Bioinfo-Master-Thesis/Output/combinedCN2.txt", sep = "\t", header = TRUE)
dataCopyNumber3 = read.csv("/Users/Yibing/Bioinfo-Master-Thesis/Output/combinedCN3.txt", sep = "\t", header = TRUE)

dataCopyNumber0$position = (dataCopyNumber0$end - dataCopyNumber0$start)/2 + dataCopyNumber0$start
dataCopyNumber1$position = (dataCopyNumber1$end - dataCopyNumber1$start)/2 + dataCopyNumber1$start
dataCopyNumber2$position = (dataCopyNumber2$end - dataCopyNumber2$start)/2 + dataCopyNumber2$start
dataCopyNumber3$position = (dataCopyNumber3$end - dataCopyNumber3$start)/2 + dataCopyNumber3$start

#cn=0  ????
dataCopyNumber0 = dataCopyNumber0[c(7, 9)]
logRCopyNumber0 = melt(dataCopyNumber0, "position")
ggplot(logRCopyNumber0, aes(value)) + geom_density() + labs(x = "logR_cn0", y = "Density", color = NULL) 
poisson_logR_cn0 = fitdist(logRCopyNumber0$value, 'pois', method = 'mme')   #?????methods
poisson_logR_cn0
#cn=1
dataCopyNumber1 = dataCopyNumber1[c(7, 9)]
logRCopyNumber1 = melt(dataCopyNumber1, "position")
ggplot(logRCopyNumber1, aes(value)) + geom_density() + labs(x = "logR_cn1", y = "Density", color = NULL) 
poisson_logR_cn1 = fitdist(logRCopyNumber1$value, 'pois', method = 'mme')   #?????methods
poisson_logR_cn1
#cn=2
dataCopyNumber2 = dataCopyNumber2[c(7, 9)]
logRCopyNumber2 = melt(dataCopyNumber2, "position")
ggplot(logRCopyNumber2, aes(value)) + geom_density() + labs(x = "logR_cn2", y = "Density", color = NULL) 
poisson_logR_cn2 = fitdist(logRCopyNumber2$value, 'pois', method = 'mme')   #?????methods
poisson_logR_cn2
#cn=3
dataCopyNumber3 = dataCopyNumber3[c(7, 9)]
logRCopyNumber3 = melt(dataCopyNumber3, "position")
ggplot(logRCopyNumber3, aes(value)) + geom_density() + labs(x = "logR_cn3", y = "Density", color = NULL) 
poisson_logR_cn3 = fitdist(logRCopyNumber3$value, 'pois', method = 'mme')   #?????methods
poisson_logR_cn3

#CORRECTED!!! separate copy number(0, 1, 2, 3) vs logR (without #reads/total reads)
dataCopyNumber0_corrected = read.csv("/Users/Yibing/Bioinfo-Master-Thesis/Output/combinedCN0_corrected.txt", sep = "\t", header = TRUE)
dataCopyNumber1_corrected = read.csv("/Users/Yibing/Bioinfo-Master-Thesis/Output/combinedCN1_corrected.txt", sep = "\t", header = TRUE)
dataCopyNumber2_corrected = read.csv("/Users/Yibing/Bioinfo-Master-Thesis/Output/combinedCN2_corrected.txt", sep = "\t", header = TRUE)
dataCopyNumber3_corrected = read.csv("/Users/Yibing/Bioinfo-Master-Thesis/Output/combinedCN3_corrected.txt", sep = "\t", header = TRUE)

dataCopyNumber0_corrected$position = (dataCopyNumber0_corrected$end - dataCopyNumber0_corrected$start)/2 + dataCopyNumber0_corrected$start
dataCopyNumber1_corrected$position = (dataCopyNumber1_corrected$end - dataCopyNumber1_corrected$start)/2 + dataCopyNumber1_corrected$start
dataCopyNumber2_corrected$position = (dataCopyNumber2_corrected$end - dataCopyNumber2_corrected$start)/2 + dataCopyNumber2_corrected$start
dataCopyNumber3_corrected$position = (dataCopyNumber3_corrected$end - dataCopyNumber3_corrected$start)/2 + dataCopyNumber3_corrected$start
#cn=0  ????
dataCopyNumber0_corrected = dataCopyNumber0_corrected[c(9, 10)]
logRCopyNumber0_corrected = melt(dataCopyNumber0_corrected, "position")
plot_cn0 = ggplot(logRCopyNumber0_corrected, aes(value)) + geom_density() + labs(x = "logR_corrected_cn0", y = "Density", color = NULL) 
poisson_logR_cn0 = fitdist(logRCopyNumber0_corrected$value, 'pois', method = 'mme')   #?????methods
poisson_logR_cn0

#cn=1
dataCopyNumber1_corrected = dataCopyNumber1_corrected[c(9, 10)]
logRCopyNumber1_corrected = melt(dataCopyNumber1_corrected, "position")
plot_cn1 = ggplot(logRCopyNumber1_corrected, aes(value)) + geom_density() + labs(x = "logR_corrected_cn1", y = "Density", color = NULL) 
logRCopyNumber1_corrected$value = logRCopyNumber1_corrected$value + 1
poisson_logR_cn1 = fitdist(logRCopyNumber1_corrected$value, "pois", method = 'mme')   #?????methods
poisson_logR_cn1 #0.9999388

#cn=2
dataCopyNumber2_corrected = dataCopyNumber2_corrected[c(9, 10)]
logRCopyNumber2_corrected = melt(dataCopyNumber2_corrected, "position")
plot_cn2 = ggplot(logRCopyNumber2_corrected, aes(value)) + geom_density() + labs(x = "logR_corrected_cn2", y = "Density", color = NULL) 
logRCopyNumber2_corrected$value = logRCopyNumber2_corrected$value + 2
poisson_logR_cn2 = fitdist(logRCopyNumber2_corrected$value, "pois",  method = 'mme')  #?????methods
poisson_logR_cn2 #1.999995


#cn=3
dataCopyNumber3_corrected = dataCopyNumber3_corrected[c(9, 10)]
logRCopyNumber3_corrected = melt(dataCopyNumber3_corrected, "position")
plot_cn3 = ggplot(logRCopyNumber3_corrected, aes(value)) + geom_density() + labs(x = "logR_corrected_cn3", y = "Density", color = NULL) 
dataCopyNumber3_corrected$value = dataCopyNumber3_corrected$value + 3
poisson_logR_cn3 = fitdist(logRCopyNumber3_corrected$value, 'pois', method = 'mle')   #?????methods
poisson_logR_cn3 #0.7042147

grid.arrange(plot_cn0, plot_cn1, plot_cn2, plot_cn3, ncol = 2, nrow = 2)
###### fitdist method ????? no maximum likelihood?#################
#############################################################################################
#Poisson
install.packages("fitdistrplus")
library(fitdistrplus)

###################
data_47 = read.csv("/Users/Yibing/Bioinfo-Master-Thesis/Output/TVEMB47.txt", sep = "\t", header = TRUE)
data_47$position = (data_47$End - data_47$Start)/2 + data_47$Start
columnNumber_cn = c(seq(7,ncol(data_47)-2,3), ncol(data_47))
data_47_cn = data_47[columnNumber_cn]
cn = (melt(data_47_cn, "position"))
cn = cn[order(cn$position),]
plotData = data.frame(filter(cn, (cn$position==unique(cn$position)[1]))) #bin #1
ggplot(plotData, aes(value)) + geom_density() + labs(x = "bin", y = "Density", color = NULL) 


poisson_cn = fitdist(plotData$value, 'pois', method = 'mle')
print(poisson_cn)

