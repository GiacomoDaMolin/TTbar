#!/usr/bin/bash

# read the input parameters
# $1: the path of the dataset
# $2: the path of the output folder 
# $3: the cross-section of the dataset
# $4: IntLumi of the dataset

# sourde the CMSSW environment
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/g/gdamolin/CMSSW_12_4_1_patch1/src && cmsenv && cd /eos/user/g/gdamolin/TT2bbemu || exit
voms-proxy-init --voms cms --valid 300:00
cp /tmp/x509up_u151129 /afs/cern.ch/user/g/gdamolin/private
# exit if cd fails
datafiles=$(dasgoclient -query="file dataset=$1")
if [ ! -d "$2$1" ]; then
    mkdir -p $2$1 || exit 
fi
for file in $datafiles; do
    #create the name of the output file: strip the path
    filename=$(echo $file | sed 's|\(^.*/\)\([a-z,A-Z,0-9,-]*\).root$|\2|')
    ofilename=$2$1/$filename"_MA".root
    echo "./MC_SignalAnalysis root://cms-xrd-global.cern.ch//$file $ofilename $3 $4" 
    condor_submit SubmitCondorMC.sub "root://cms-xrd-global.cern.ch//$file $ofilename $3 $4" 
done
