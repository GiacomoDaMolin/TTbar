#!/usr/bin/bash

echo "start"
usage() { echo "Usage: $0 [-f <submitfile>] [-e <executable> ] [-d <dataset>] [-o <outpath>] [-x xsec] [-l <lumi>] [-p <user_proxy>]" 1>&2; exit 1; }
while getopts "f:j:d:o:x:l:p:" opt; do
    case "$opt" in
        f) SUBFILE=$OPTARG
            ;;
        j) JOBFILE=$OPTARG
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
echo "Using dasgoclient to list files and create submission jobs"
datafiles=`dasgoclient -query="file dataset=${DATASET}"`
jobdescription=${SUBFILE}
if [[ -f $jobdescription ]]; then
    rm ${jobdescription}
fi
for file in ${datafiles[@]}; do
    echo "-e ${EXE} -d root://cms-xrd-global.cern.ch//${file} -o ${OUTPATH} -x ${XSEC} -l ${LUMI} -p ${X509_USER_PROXY}" >> ${jobdescription}
done
echo "Use ${jobdescription} with your condor file"
