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
    #create .txt and write arguments on .txt
    echo "$1, root://cms-xrd-global.cern.ch/${filein}, ${ofilename}, $4, $5, $6," > myfile.txt
    # condor submit, .sub file must read myfile.txt
    "condor_submit /afs/cern.ch/user/g/gdamolin/Johan/TTbar/TJobforCondor.sub"
    #sleep 2 seconds
    sleep 2
    #delete .txt
    rm myfile.txt
    echo "##################"
done
