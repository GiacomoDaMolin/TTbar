# Workflow
Selections & preliminar histograms are done in:
Mixed_Analysis.cpp for MC
Data_Analysis.cpp for Data.

These .cpp file have their histograms defined globally in Histodef.cpp and call some (long but pretty easy) functions of Auxiliary.cpp
The logical structure is the foolowing
-the histograms are created globally, the main is called and ->sumw2() is applied to all histograms 
-then I set the input tree from the NanoAOD ready to be red, as well as preparing the corrections functions.
-I loop in the events, perform the selections, compute weights and corrections and finally fill them in Th1Ds.
-After the loop I just cout some numbers and Write the histo in the outfile (which will be in Eos).

To continue with the analysis, one needs to hadd the files of each sample, divide their yield by the SumOfMCweights and then plot them as he/she pleases.
I have a script for that, I can pass it to you if you want. However the names of files and directories are hard-coded.

How do I launc the Analysis?
I use createdagROOT.py (see the ReadMe of DY->MUMU) that creates a submit directory.
From there I just run condor_submit_dag for each directory (created from the .json).
If some jobs fail, one can resubmit with Python_Analysis/rescue.py, which will send in condor only the jobs that failed. (see DY->MUMU readME).



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
