#!/usr/bin/bash

# read the input parameters
# $1: the executable (for example: ./Data_Analysis)
# $2: the path of the dataset
# $3: the path of the output folder 
# $4: the cross-section of the dataset
# $5: IntLumi of the dataset
# $6: the signal/bkgd flag
# $7: proxy path

# exit if cd fails
export X509_USER_PROXY=/afs/cern.ch/user/g/gdamolin/private/x509up_u151129
#voms-proxy-info -all
voms-proxy-info -all -file /afs/cern.ch/user/g/gdamolin/private/x509up_u151129

echo "Before"
echo "$2"
echo "This is my argument"
echo "$1 $3 $4 $5 $6"
echo "These are the rest"

datafiles=$(dasgoclient -query="file dataset=root://cms-xrd-global.cern.ch///store/mc/RunIISummer20UL18NanoAODv2$2")

#datafiles=$(dasgoclient -query="file dataset=$2")

echo $datafiles

for file in $datafiles; do
    #create the name of the output file: strip the path
    filename=$(echo $file | sed 's|\(^.*/\)\([a-z,A-Z,0-9,-]*\).root$|\2|')
    echo $filename
    ofilename=$3/$filename"_MA".root
    echo "$1 root://cms-xrd-global.cern.ch//$file $ofilename $4 $5 $6" 
    "$1 root://cms-xrd-global.cern.ch//$file $ofilename $4 $5 $6" 
done
