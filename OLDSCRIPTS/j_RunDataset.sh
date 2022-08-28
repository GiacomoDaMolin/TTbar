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
voms-proxy-info -all
#voms-proxy-info -all -file /afs/cern.ch/user/g/gdamolin/private/x509up_u151129

datafiles=$(dasgoclient -query="file dataset=/DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1/NANOAODSIM")

echo "$datafiles"

for file in $datafiles; do
    #create the name of the output file: strip the path
    filename=$(echo $file | sed 's|\(^.*/\)\([a-z,A-Z,0-9,-]*\).root$|\2|')
    ofilename="/afs/cern.ch/user/g/gdamolin/Johan/$filename_MA.root"
    echo "./Mixed_Analysis root://cms-xrd-global.cern.ch//$file $ofilename 951.5*1000 59.83 False"
    #root 
    "./Mixed_Analysis root://cms-xrd-global.cern.ch//$file $ofilename 951.5*1000 59.83 False"
    #.q 
done
