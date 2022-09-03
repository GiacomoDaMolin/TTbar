# How to send jobs in Condor

To create the description of the jobs needed to run over a dataset the following command can be used

```
sh RunDataset.sh -e myanalysis.py \
    -d dataset_name \
    -o output_dir
    -x xsec \
    -l integ_luminosity \
    -s is_signal
    -p full_path_to_proxy \
    -c full_path_to_cmssw
```
NOTE: use RunDataset.sh if you want to run on MC, RunActualData.sh if you wish to run on Data.

The output is a `jobdescription.txt` file which will list all the arguments to be used to process single files.
This will be called when you launch the job on condor.
Before submitting to condor, test locally that a single file will be processed as desired by doing

```
sh RunDataset.sh `head -n 1 jobdescription.txt`
```

If all looks good you can submit to condor by doing

```
condor_submit jobdescription.sub
```

LAWS OF CONDOR:
1) Thou cannot put your logs/executable nor output files in eos, Condor will refuse to read it. You can however save them in afs and then move them in eos in the script.
2) Thou cannot use mkdir while using condor.
3) Thou need to create the c++ executable with THE SAME CMSSW CMSENV version of the one you put in full_path_to_cmssw
4) I will share soon a file with all the commands for each dataset, possibly we can make a script from it
5) Thou will help me make me a Postscript that "hadd"s all the output files in a single file because doing that manually is a pain



# Data 

## EGamma
### A-C are v9-v1
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FEGamma%2FRun2018A-UL2018_MiniAODv2_NanoAODv9-v1%2FNANOAOD&instance=prod/global
### D is v9-v3
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FEGamma%2FRun2018D-UL2018_MiniAODv2_NanoAODv9-v3%2FNANOAOD&instance=prod/global

### Luminosity

```
brilcalc lumi -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json --hltpath "HLT_Ele32_WPTight_Gsf_v*" -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt
```

#Summary:  
+---------------------------+-------+------+--------+-------------------+------------------+  
| hltpath                   | nfill | nrun | ncms   | totdelivered(/fb) | totrecorded(/fb) |  
+---------------------------+-------+------+--------+-------------------+------------------+  
| HLT_Ele32_WPTight_Gsf_v13 | 23    | 50   | 22670  | 5.449751450       | 5.291407965      |  
| HLT_Ele32_WPTight_Gsf_v14 | 10    | 26   | 13311  | 3.796985905       | 3.659410871      |  
| HLT_Ele32_WPTight_Gsf_v15 | 164   | 402  | 198129 | 52.965171486      | 50.877060670     |  
+---------------------------+-------+------+--------+-------------------+------------------+  
#Sum delivered : 62.211908842  
#Sum recorded : 59.827879506  
#Check JSON:  
#(run,ls) in json but not in results: [(325172, 477), (325172, 478), (325172, 479), (325172, 480), (325172, 481), (325172, 482), (325172, 483), (325172, 484), (325172, 485)]  

## SingleMuon
### A-C are v9-v2
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FSingleMuon%2FRun2018A-UL2018_MiniAODv2_NanoAODv9-v2%2FNANOAOD&instance=prod/global
### D is v9-v1
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FSingleMuon%2FRun2018D-UL2018_MiniAODv2_NanoAODv9-v1%2FNANOAOD&instance=prod/global

### Luminosity:
```
brilcalc lumi -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json --hltpath "HLT_IsoMu24_v*" -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt
```
#Summary:  
+-----------------+-------+------+--------+-------------------+------------------+  
| hltpath         | nfill | nrun | ncms   | totdelivered(/fb) | totrecorded(/fb) |  
+-----------------+-------+------+--------+-------------------+------------------+  
| HLT_IsoMu24_v11 | 23    | 50   | 22670  | 5.449751450       | 5.291407965      |  
| HLT_IsoMu24_v12 | 10    | 26   | 13287  | 3.788621855       | 3.651245839      |  
| HLT_IsoMu24_v13 | 164   | 402  | 198129 | 52.965171486      | 50.877060670     |  
+-----------------+-------+------+--------+-------------------+------------------+  
#Sum delivered : 62.203544792  
#Sum recorded : 59.819714474  
#Check JSON:  
#(run,ls) in json but not in results: [(325172, 477), (325172, 478), (325172, 479), (325172, 480), (325172, 481), (325172, 482), (325172, 483), (325172, 484), (325172, 485)]  
# MC

## Signal
### Fully Leptonic Dataset 
GenCrossSection=72.1pb
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FTTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global

### (unused) Semi-Leptonic Dataset 
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FTTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global


## Backgrounds 
### Drell-Yann (1Jet) Background production.
GenCrossSection=951.5 pb (or 955.8)
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FDYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global
```
./run_MC.sh ./Mixed_Analysis /DYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1/NANOAODSIM ./output/MC 951.5*1000 59.83 False
```

### DY 2 Jets
GenCrossSection=360.4 pb 
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FDYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global
```
./run_MC.sh ./Mixed_Analysis /DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1/NANOAODSIM ./output/MC/ 360.4*1000 59.83 False
```

### DY Mass Range
GenCrossSection=15810 pb
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FDYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global
```
./run_MC.sh ./Mixed_Analysis /DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1/NANOAODSIM ./output/MC/ 15810*1000 59.83 False
```
### W 0 Jet
GenCrossSection=53330 pb, #Events: 177728297
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FWJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global
```
./run_MC.sh ./Mixed_Analysis /WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1/NANOAODSIM ./output/MC 53330*1000 59.83 False
```

### W 1 Jet
GenCrossSection=8875.0 pb, #Events: 190471282
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FWJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global
```
./run_MC.sh ./Mixed_Analysis /WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1/NANOAODSIM ./output/MC 8875.0*1000 59.83 False
```
### W 2 Jet
GenCrossSection=3338.0 pb
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FWJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global
```
./run_MC.sh ./Mixed_Analysis /WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1/NANOAODSIM ./output/MC 3338 59.83*1000 False
```

### SingleTop
GenCrossSection=32.45 pb
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV_PDFWeights-powheg-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global
```
./run_MC.sh ./Mixed_Analysis /ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV_PDFWeights-powheg-pythia8/RunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1/NANOAODSIM ./output/MC 32.45*1000 59.83 False
```
### SingleAntiTop
GenCrossSection=32.51 pb
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV_PDFWeights-powheg-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global
```
./run_MC.sh ./Mixed_Analysis /ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV_PDFWeights-powheg-pythia8/RunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1/NANOAODSIM ./output/MC 32.51*1000 59.83 False
```
