install.packages("dplyr")
library(readr)
library(dplyr)
library(reshape2)
library(ggplot2)

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

  pdf(plotDataName_logR)
  plot = ggplot(plotDataCombined_logR, aes(position, value, group=variable, color = variable)) +
    geom_point(size = 0.2) + 
    facet_wrap(~Chr, scales = "free_x", ncol=3)+
    labs(x = "Position", y = "logR", color = NULL) +
    theme(legend.position = "bottom")
  print(plot)
  dev.off()
  
  png(plotDataName_cn)
  plot = ggplot(plotDataCombined_cn, aes(position, value, group=variable, color = variable)) +
    geom_point(size = 0.2) + 
    facet_wrap(~Chr, scales = "free_x", ncol=3)+
    labs(x = "Position", y = "cn", color = NULL) +
    theme(legend.position = "bottom")
  print(plot)
  dev.off()
  
}

#### order chromosome!!!!
#### x y axis scale!!!!
####graph region enlarge!!!!
############################################################################

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
ggplot(plotData, aes(value)) + geom_density() + labs(x = "cn", y = "Density", color = NULL) 


poisson_cn = fitdist(plotData$value, 'pois', method = 'mme')
print(poisson_cn)
###### fitdist method ????? no maximum likelihood?
