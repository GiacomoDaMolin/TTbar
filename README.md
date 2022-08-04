# Data 

## EGamma
### A-C are v9-v1
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FEGamma%2FRun2018A-UL2018_MiniAODv2_NanoAODv9-v1%2FNANOAOD&instance=prod/global
### D is v9-v3
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FEGamma%2FRun2018D-UL2018_MiniAODv2_NanoAODv9-v3%2FNANOAOD&instance=prod/global

### Luminosity

```brilcalc lumi -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json --hltpath "HLT_Ele32_WPTight_Gsf_v*" -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt```

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
```brilcalc lumi -u /fb --normtag /cvmfs/cms-bril.cern.ch/cms-lumi-pog/Normtags/normtag_PHYSICS.json --hltpath "HLT_IsoMu24_v*" -i /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions18/13TeV/Legacy_2018/Cert_314472-325175_13TeV_Legacy2018_Collisions18_JSON.txt```
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

### DY 2 Jets
GenCrossSection=360.4 pb 
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FDYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global

### DY Mass Range
GenCrossSection=15810 pb
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FDYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global
### W 0 Jet
GenCrossSection=53330 pb, #Events: 177728297
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FWJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global

### W 1 Jet
GenCrossSection=8875.0 pb, #Events: 190471282
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FWJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global
### W 2 Jet
GenCrossSection=3338.0 pb
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FWJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global

### SingleTop
GenCrossSection=32.45 pb
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV_PDFWeights-powheg-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global
### SingleAntiTop
GenCrossSection=32.51 pb
https://cmsweb.cern.ch/das/request?input=dataset%3D%2FST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV_PDFWeights-powheg-pythia8%2FRunIISummer20UL18NanoAODv2-106X_upgrade2018_realistic_v15_L1v1-v1%2FNANOAODSIM&instance=prod/global