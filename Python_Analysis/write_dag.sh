#!/usr/bin/bash

echo "start"
usage() { echo "Usage: $0 [-f <submitfile>] [-j <.dag jobfile> ] [-d <dataset>] [-o <outpath>] [-x xsec] [-l <lumi>] [-p <user_proxy>]" 1>&2; exit 1; }
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
if [[ -f ${JOBFILE} ]]; then
    rm ${JOBFILE}
fi
for file in ${datafiles[@]}; do
    JOBID=(${file//\// })
    JOBID=${JOBID[-1]}
    echo "JOB ${JOBID} ${SUBFILE}" >> ${JOBFILE}
    echo "VARS ${JOBID} INFILE=\"root://cms-xrd-global.cern.ch//${file}\" OUTFILE=\"${OUTPATH}\" XS=\"${XSEC}\" LUMI=\"${LUMI}\" PROXY=\"${X509_USER_PROXY}\"" >> ${JOBFILE}
done
echo "execute jobs with 'condor_submit_dag -config dagman.config ${JOBFILE}'"
echo "or just:"
echo "execute jobs with 'condor_submit_dag ${JOBFILE}'"
echo "remember to put ${SUBFILE} in the same dir"
