library(readr)

getwd()
setwd("/Users/Yibing/Desktop/Leuven/Thesis")

#raw data
rawDataHeader = c("chr", "start", "end", "binSize", "length", "readNumber", "???")
TVEMB1_1_250K_101_sorted = read.csv("/Users/Yibing/Desktop/Leuven/Thesis/data/raw_counts/TVEMB1_1.250K.101.sorted.count-t", sep = "\t", header = FALSE)
TVEMB1_2.250K.101.sorted = read.csv("/Users/Yibing/Desktop/Leuven/Thesis/data/raw_counts/TVEMB1_2.250K.101.sorted.count-t", sep = "\t", header = FALSE)
colnames(TVEMB1_1_250K_101_sorted) = rawDataHeader
colnames(TVEMB1_2.250K.101.sorted) = rawDataHeader
View(TVEMB1_1_250K_101_sorted)
View(TVEMB1_2.250K.101.sorted)

#logR cn per bin
TVEMB1_1_250K_101_sorted_count_gamma_15_gc_corrected_segments_copynumber = read.csv("Desktop/Leuven/Thesis/data_for_yibing/logr_cn_per_bin/TVEMB1_1.250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.txt", sep = "\t")
View(TVEMB1_1_250K_101_sorted_count_gamma_15_gc_corrected_segments_copynumber)

#segmented cell
TVEMB1_1_250K_101_sorted_count_gamma_15_gc_corrected_segments_copynumber_breakpoints = read.csv("Desktop/Leuven/Thesis/data_for_yibing/segmented/TVEMB1_1.250K.101.sorted.count.gamma_15.gc_corrected.segments.copynumber.breakpoints.txt", sep = "\t")
View(TVEMB1_1_250K_101_sorted_count_gamma_15_gc_corrected_segments_copynumber_breakpoints)


