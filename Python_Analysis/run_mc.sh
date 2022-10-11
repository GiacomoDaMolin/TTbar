#!/usr/bin/bash

usage() { echo "Usage: $0 [-e <executable> ] [-d <dataset>] [-o <outpath>] [-x xsec] [-l <lumi>] [-w <sum_w>] [-p <user_proxy>]" 1>&2; exit 1; }
while getopts "i:o:x:l:w:p:" opt; do
    case "$opt" in
        i) INFILE=$OPTARG
            ;;
        o) OUTPATH=$OPTARG
            ;;
        x) XSEC=$OPTARG
            ;;
        l) LUMI=$OPTARG
            ;;
        w) SUM_W=$OPTARG
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
${EXE} -i $filename -o $ofilename -x $XSEC -l $LUMI -s $SUM_W -m  || { 
    echo "${EXE} failed with file ${filename}, removing intermediate file" 1>&2; 
    if [[ -f $ofilename ]]; then
    rm "$ofilename"
    fi
    exit 1
    }
mv $ofilename ${OUTPATH}/${filestring}_MA.root
