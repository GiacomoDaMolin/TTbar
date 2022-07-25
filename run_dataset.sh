#!/usr/bin/bash

# read the input parameters
# $1: the path of the dataset
# $2: the path of the output folder (no trailing "/")
# $3: the cross-section of the dataset
# $4: IntLumi of the dataset
# $5: the signal/bkgd flag 

# sourde the CMSSW environment
source /cvmfs/cms.cern.ch/cmsset_default.sh
voms-proxy-init --rfc --voms cms -valid 192:00
# exit if cd fails
cd /eos/user/j/jowulff/TTbar/CMSSW_12_4_0/src && cmsenv && cd /eos/user/j/jowulff/TTbar || exit
datafiles=$(dasgoclient -query="file dataset=$1")
for file in $datafiles; do
    #create the name of the output file: strip the path
    filename=$(echo $file | sed 's|\(^.*/\)\([a-z,A-Z,0-9,-]*\).root$|\2|')
    ofilename=$2/$filename"_MA".root
    echo "./Mixed_Analysis root://cms-xrd-global.cern.ch//$file $ofilename $3 $4 $5"
    ./Mixed_Analysis root://cms-xrd-global.cern.ch//$file $ofilename $3 $4 $5
done
