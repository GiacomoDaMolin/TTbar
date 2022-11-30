
#!/usr/bin/bash

# parse the command line arguments
# -e the executable (for example: ./Data_Analysis)
# -d the dataset name as in DAS
# -o the path of the output folder 
# -x the cross-section of the dataset
# -l IntLumi of the dataset
# -p proxy path

echo "start"
usage() { echo "Usage: $0 [-e <executable> ] [-p <dataset>] [-o <outpath>] [-p <user_proxy>] [-f <first_data>]" 1>&2; exit 1; }
while getopts "e:d:o:p:f:" opt; do
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
        f) FIRST_DATA=$OPTARG
        ;;
        *)
        usage
        ;;
    esac
done
echo "Using dasgoclient to list files and create submission jobs"
datafiles=`dasgoclient -query="file dataset=${DATASET}"`
jobdescription="jobdescription.txt"
if [[ -f $jobdescription ]]; then
    rm ${jobdescription}
fi
for file in ${datafiles[@]}; do
    echo "-e ${EXE} -d root://cms-xrd-global.cern.ch//${file} -o ${OUTPATH} -p ${X509_USER_PROXY} -f {$FIRST_DATA}" >> ${jobdescription}
done
echo "Use ${jobdescription} with your condor file"