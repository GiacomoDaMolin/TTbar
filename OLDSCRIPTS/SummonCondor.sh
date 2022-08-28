#!/usr/bin/bash

# read the input parameters
# $1: the executable (for example: ./Data_Analysis)
# $2: the path of the dataset
# $3: the path of the output folder 
# $4: the cross-section of the dataset
# $5: IntLumi of the dataset
# $6: the signal/bkgd flag


# exit if cd fails
export X509_USER_PROXY=/afs/cern.ch/user/g/gdamolin/private/x509up_u151129
voms-proxy-info -all
#voms-proxy-info -all -file /afs/cern.ch/user/g/gdamolin/private/x509up_u151129

datafiles=$(dasgoclient -query="file dataset=$2")

for filein in $datafiles; do
    #create the name of the output file: strip the path
    filename=$(echo $filein | sed 's|\(^.*/\)\([a-z,A-Z,0-9,-]*\).root$|\2|')
    echo $filein
    ofilename=$3/$filename"_MA".root
    echo "condor_submit /afs/cern.ch/user/g/gdamolin/Johan/TTbar/JobforCondor.sub $1 root://cms-xrd-global.cern.ch/${filein} ${ofilename} $4 $5 $6"
    "condor_submit /afs/cern.ch/user/g/gdamolin/Johan/TTbar/JobforCondor.sub $1 root://cms-xrd-global.cern.ch/${filein} ${ofilename} $4 $5 $6"
    #"condor_submit /afs/cern.ch/user/g/gdamolin/Johan/TTbar/JobforCondor.sub /afs/cern.ch/user/g/gdamolin/Johan/TTbar/Mixed_Analysis root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/130000/6E2A3FD3-C3FD-5A4C-B97C-58DC2C0A0A60.root /eos/user/g/gdamolin/TT2bbemu/DY1/6E2A3FD3-C3FD-5A4C-B97C-58DC2C0A0A60_MA.root 951.5*1000 59.83 False"
    echo "##################"
done
