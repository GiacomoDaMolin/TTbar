#!/usr/bin/bash

# read the input parameters
# $1: the executable (for example: ./Data_Analysis)
# $2: the path of the dataset
# $3: the path of the output folder 
# $4: the cross-section of the dataset
# $5: IntLumi of the dataset
# $6: the signal/bkgd flag

# exit if cd fails
datafiles=$(dasgoclient -query="file dataset=$1")

for file in $datafiles; do
    #create the name of the output file: strip the path
    filename=$(echo $file | sed 's|\(^.*/\)\([a-z,A-Z,0-9,-]*\).root$|\2|')
    ofilename=$3$2/$filename"_MA".root
    echo "$1 root://cms-xrd-global.cern.ch//$file $ofilename $4 $5 $6"
    root 
    ".$1 root://cms-xrd-global.cern.ch//$file $ofilename $4 $5 $6"
    .q 
done
