#!/usr/bin/bash

# parse the command line arguments
# -e the executable (for example: ./Data_Analysis)
# -d the dataset name as in DAS
# -o the path of the output folder 
# -x the cross-section of the dataset
# -l IntLumi of the dataset
# -s the signal/bkgd flag
# -p proxy path

echo "start"
X509_USER_PROXY=/afs/cern.ch/user/g/gdamolin/private/x509up_u151129
CMSSW=/afs/cern.ch/user/g/gdamolin/CMSSW_12_4_1_patch1/src
while getopts "e:d:o:x:l:s:p:c:" opt; do
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
        s) SIGNAL=$OPTARG
            ;;
        p) X509_USER_PROXY=$OPTARG
            ;;
        c) CMSSW=$OPTARG
            ;;
    esac
done

echo "stuck in function definition?"

#define a simple function to test the inputs are available
test_input () {
  echo "$1 = $2"
  if [ -z "$2" ]
  then
      echo "$1 is empty, use option -$3 to define... quitting"
      exit -1
  fi
}

echo "test2"
test_input "Executable" "$EXE" "e"
test_input "Dataset" "$DATASET" "d"
test_input "Output path" "$OUTPATH" "o"
test_input "Cross section" "$XSEC" "x"
test_input "Integ. luminosity" "$LUMI" "l"
test_input "Signal" "$SIGNAL" "s"
test_input "Path to proxy" "$X509_USER_PROXY" "p"
test_input "CMSSW" "$CMSSW" "c"

#move file to eos
outdir="/eos/user/g/gdamolin/TT2bbemu/temp"

#now test what has been asked to analyze
#1. a single file run executable
#2. a full dataset, create a condor file
if [[ "$DATASET" == *"root://"* ]]; then

  echo "Dataset is already a single file, will run executable for it"

  filename=$DATASET
  filestring=$(echo $filename | sed 's|\(^.*/\)\([a-z,A-Z,0-9,-]*\).root$|\2|')
  ofilename=${OUTPATH}/$filestring"_MA".root

  cd $CMSSW
  eval `scram r -sh`
  export X509_USER_PROXY=$X509_USER_PROXY
  cd -

  echo "CMSSW_BASE is now set to $CMSSW_BASE"
  echo "PROXY is now set to $X509_USER_PROXY"
  echo "Executing analysis script as"
  ${EXE} $filename $ofilename ${XSEC} ${LUMI} ${SIGNAL}
  mv $ofilename ${outdir}/${filestring}_MA.root

else

  echo "Using dasgoclient to list files and create submission jobs"

  datafiles=`dasgoclient -query="file dataset=${DATASET}"`
  jobdescription="jobdescription.txt"
  rm ${jobdescription}
  for file in ${datafiles[@]}; do
      echo "-e ${EXE} -d root://cms-xrd-global.cern.ch//${file} -o ${OUTPATH} -x ${XSEC} -l ${LUMI} -s ${SIGNAL} -p ${X509_USER_PROXY} -c $CMSSW" >> ${jobdescription} 
  done
  echo "Use ${jobdescription} with your condor file"

fi


