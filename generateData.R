args = commandArgs()

#normal distribution parameters
mean_cn1 = -0.8269323
mean_cn2 = 0.01823383
mean_cn3 = 0.5264469
sd_cn1 = 0.5603306
sd_cn2 = 0.41483619
sd_cn3 = 0.4293003

Single_CNV = read.table("single_cnv.txt", sep="\t", header = TRUE)

row_number = as.numeric(args[3])

colnames(Single_CNV) = c("total_number_of_samples", "length", "cn2_number_of_samples", "cn2_proportion", "cn1_number_of_samples", "cn1_proportion", "cnv_size", "size_proportion", "cn1_start", "cn1_end",  "cn3_number_of_samples", "cn3_proportion", "cnv_size", "size_proportion", "cn3_start", "cn3_end")
length = 1000
cn1_start_position = as.numeric(Single_CNV$cn1_start[row_number])
cn3_start_position = as.numeric(Single_CNV$cn3_start[row_number])
cnv_size = as.numeric(Single_CNV$cnv_size[row_number])
total_number = as.numeric(Single_CNV$total_number_of_samples[row_number])
number_no_cnv = as.numeric(Single_CNV$cn2_number_of_samples[row_number])
number_cnv_cn1 = as.numeric(Single_CNV$cn1_number_of_samples[row_number])
number_cnv_cn3 = as.numeric(Single_CNV$cn3_number_of_samples[row_number])

for (n in 1:10){
  stimulate_logR = data.frame("position" = 1:length)
  #no cnv stimulated data
  for (i in 1:number_no_cnv){
    column_name = paste("sample", i, sep = "_")
    stimulate_sample = as.data.frame(rnorm(length, mean = mean_cn2, sd = sd_cn2))
    colnames(stimulate_sample) = column_name
    stimulate_logR = cbind(stimulate_logR, stimulate_sample)
  }
  
  #cn1 stimulated data
  if (number_cnv_cn1 != 0){
    for (i in 1:number_cnv_cn1){
      column_name = paste("sample", number_no_cnv+i, sep = "_")
      cnv_cn = 1
      if (cn1_start_position!=1){
        stimulate_segment_before_cnv = data.frame("position" = 1:(cn1_start_position-1), "reads" = rnorm((cn1_start_position-1), mean = mean_cn2, sd = sd_cn2))
        stimulate_segment_cnv = data.frame("position" = cn1_start_position:(cn1_start_position+cnv_size-1), "reads" = rnorm((cnv_size), mean = mean_cn1, sd = sd_cn1))
        stimulate_segment_after_cnv = data.frame("position" = (cn1_start_position+cnv_size):length, "reads" = rnorm((length-cn1_start_position-cnv_size+1), mean = mean_cn2, sd = sd_cn2))
        stimulate_sample = rbind(stimulate_segment_before_cnv, stimulate_segment_cnv, stimulate_segment_after_cnv)[2]
        colnames(stimulate_sample) = column_name
        stimulate_logR = cbind(stimulate_logR, stimulate_sample)
      }
      else{
        stimulate_segment_cnv = data.frame("position" = cn1_start_position:(cn1_start_position+cnv_size-1), "reads" = rnorm((cnv_size), mean = mean_cn1, sd = sd_cn1))
        stimulate_sample = rbind(stimulate_segment_cnv)[2]
        colnames(stimulate_sample) = column_name
        stimulate_logR = cbind(stimulate_logR, stimulate_sample)
      }
    }
  }
  #cn3 stimulated data
  if (number_cnv_cn3 != 0){
    for (i in 1:number_cnv_cn3){
      column_name = paste("sample", number_no_cnv+number_cnv_cn1+i, sep = "_")
      cnv_cn = 3
      if (cn3_start_position!=1){
        stimulate_segment_before_cnv = data.frame("position" = 1:(cn3_start_position-1), "reads" = rnorm((cn3_start_position-1), mean = mean_cn2, sd = sd_cn2))
        stimulate_segment_cnv = data.frame("position" = cn3_start_position:(cn3_start_position+cnv_size-1), "reads" = rnorm((cnv_size), mean = mean_cn3, sd = sd_cn3))
        stimulate_segment_after_cnv = data.frame("position" = (cn3_start_position+cnv_size):length, "reads" = rnorm((length-cn3_start_position-cnv_size+1), mean = mean_cn2, sd = sd_cn2))
        stimulate_sample = rbind(stimulate_segment_before_cnv, stimulate_segment_cnv, stimulate_segment_after_cnv)[2]
        colnames(stimulate_sample) = column_name
        stimulate_logR = cbind(stimulate_logR, stimulate_sample)
      }
      else{
        stimulate_segment_cnv = data.frame("position" = cn3_start_position:(cn3_start_position+cnv_size-1), "reads" = rnorm((cnv_size), mean = mean_cn3, sd = sd_cn3))
        stimulate_sample = rbind(stimulate_segment_cnv)[2]
        colnames(stimulate_sample) = column_name
        stimulate_logR = cbind(stimulate_logR, stimulate_sample)
      }
    }
  }
  writeFileName = paste(paste(paste(paste(paste(paste(paste(paste("/Users/Yibing/Bioinfo-Master-Thesis/Thesis/generateData/SD_",total_number,sep = ""), number_no_cnv,sep = "_"), number_cnv_cn1, sep = "_"), number_cnv_cn3, sep = "_"), cnv_size, sep = "_"), "run", sep = "_"), n, sep = ""), "csv", sep = ".")
  write.csv(stimulate_logR, file = writeFileName)
}
