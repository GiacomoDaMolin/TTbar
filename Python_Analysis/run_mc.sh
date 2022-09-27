#!/usr/bin/bash

usage() { echo "Usage: $0 [-e <executable> ] [-p <dataset>] [-o <outpath>] [-x xsec] [-l <lumi>] [-p <user_proxy>]" 1>&2; exit 1; }
while getopts "e:d:o:x:l:p:f:" opt; do
    case "$opt" in
        e) EXE=$OPTARG
            ;;
        d) DATASET=$OPTARG
            ;;
        o) OUTPATH=$OPTARG
            ;;
        x) XSEC=$OPTARG
            ;;
        l) LUMI=$OPTARG
            ;;
        p) X509_USER_PROXY=$OPTARG
            ;;
        *)
        usage
        ;;
    esac
done
outdir="/afs/cern.ch/user/j/jowulff/Condor/TTbar/MC"
filename=$DATASET
filestring=$(echo $filename | sed 's|\(^.*/\)\([a-z,A-Z,0-9,-]*\).root$|\2|')
ofilename=${outdir}/$filestring"_MA".root

export X509_USER_PROXY=$X509_USER_PROXY

echo "PROXY is now set to $X509_USER_PROXY"
${EXE} -i $filename -o $ofilename -x $XSEC -l $LUMI -m 
mv $ofilename ${OUTPATH}/${filestring}_MA.root
