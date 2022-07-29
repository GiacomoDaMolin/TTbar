#!/usr/bin/bash
usage() { echo "Usage: $0 [-p <path to dataset>] [-o <path to output folder>]" 1>&2; exit 1; }

while getopts "p:o:" i; do
    case "${i}" in
        p)
            p=${OPTARG}
            ;;
        o)
            o=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done
if [ -z "${p}" ] || [ -z "${o}" ]; then
   usage
fi
# sourde the CMSSW environment
source /cvmfs/cms.cern.ch/cmsset_default.sh
voms-proxy-init --rfc --voms cms -valid 192:00
# exit if cd fails
cd /eos/user/j/jowulff/TTbar/CMSSW_12_4_0/src && cmsenv && cd /eos/user/j/jowulff/TTbar || exit
datafiles=$(dasgoclient -query="file dataset=$p")
# check if the output folder exists if not create it
if [ ! -d "$o$p" ]; then
    mkdir -p $o$p
fi
for file in $datafiles; do
    #create the name of the output file: strip the path
    filename=$(echo $file | sed 's|\(^.*/\)\(.*$\)|\2|')
    ofilename=$o$p/$filename
    echo "./Data_Analysis root://cms-xrd-global.cern.ch//$file $ofilename"
    ./Data_Analysis root://cms-xrd-global.cern.ch//$file $ofilename
done
