#!/usr/bin/bash

usage() { echo "Usage: $0 [-e <executable> ] [-d <dataset>] [-o <outpath>] [-p <user_proxy>]" 1>&2; exit 1; }
while getopts "i:o:f:p:" opt; do
    case "$opt" in
        i) INFILE=$OPTARG
            ;;
        o) OUTPATH=$OPTARG
            ;;
        f) FIRST_DATA=$OPTARG
            ;;
        p) X509_USER_PROXY=$OPTARG
            ;;
        *)
        usage
        ;;
    esac
done
outdir="/afs/cern.ch/user/j/jowulff/Condor/TTbar/MC"
filename=$INFILE
filestring=$(echo $filename | sed 's|\(^.*/\)\([a-z,A-Z,0-9,-]*\).root$|\2|')
ofilename=${outdir}/$filestring"_MA".root
EXE="/afs/cern.ch/user/j/jowulff/Condor/TTbar/PythonAnalysis/MixedAnalysis.py"

export X509_USER_PROXY=$X509_USER_PROXY

echo "PROXY is now set to $X509_USER_PROXY"
if [ "$FIRST_DATA" == true ]; then 
${EXE} -i $filename -o $ofilename --first_data || { 
    echo "${EXE} failed with file ${filename}, removing intermediate file" 1>&2; 
    if [[ -f $ofilename ]]; then
    rm "$ofilename"
    fi
    exit 1
    }
mv $ofilename ${OUTPATH}/${filestring}_MA.root
else
${EXE} -i $filename -o $ofilename || { 
    echo "${EXE} failed with file ${filename}, removing intermediate file" 1>&2; 
    if [[ -f $ofilename ]]; then
    rm "$ofilename"
    fi
    exit 1
    }
mv $ofilename ${OUTPATH}/${filestring}_MA.root
fi