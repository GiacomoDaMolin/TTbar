#!/usr/bin/bash

echo "PreScript beginning"

# sourde the CMSSW environment
source /cvmfs/cms.cern.ch/cmsset_default.sh
cd /afs/cern.ch/user/g/gdamolin/CMSSW_12_4_1_patch1/src && cmsenv && cd /eos/user/g/gdamolin/TT2bbemu || exit
#voms-proxy-init --voms cms --valid 300:00
#cp /tmp/x509up_u151129 /afs/cern.ch/user/g/gdamolin/private

pwd
echo "PreScript done"
