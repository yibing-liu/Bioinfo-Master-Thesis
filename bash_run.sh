#!/bin/bash
#generate data
for row_number in {1..1}
do
    #generate data
    R < generateData.R --no-save $row_number
    for files in "/Users/Yibing/Bioinfo-Master-Thesis/Thesis/generateData/"*.csv
    do
        /Applications/MATLAB_R2019a.app/bin/matlab -nodisplay -nodesktop -r "filename='$files';run pla.m;quit"
    done
    #run segmentation&calculate copynumber
    R < scriptName.R --no-save

#output result
done
