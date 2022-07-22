## The 019*.root file is fully leptonic
GenCrossSection=72.1pb
https://cmsweb.cern.ch/das/request?input=file%3D%2Fstore%2Fmc%2FRunIISummer20UL18NanoAODv2%2FTTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8%2FNANOAODSIM%2F106X_upgrade2018_realistic_v15_L1v1-v1%2F270000%2F019426EE-3D50-1249-B266-F6DBA0AFE3B5.root&instance=prod/global


## The 01A*.root file is semileptonic!!
https://cmsweb.cern.ch/das/request?input=file%3D%2Fstore%2Fmc%2FRunIISummer20UL18NanoAODv2%2FTTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8%2FNANOAODSIM%2F106X_upgrade2018_realistic_v15_L1v1-v1%2F20000%2F01A179D4-6A3E-BD4F-BDD8-A3B3E99D44F9.root&instance=prod/global


## The 162*.root file is a Drell-Yann (1Jet) Background production.
GenCrossSection=951.5 fb (or 955.8)
https://cmsweb.cern.ch/das/request?input=file%3D%2Fstore%2Fmc%2FRunIISummer20UL18NanoAODv2%2FDYJetsToLL_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8%2FNANOAODSIM%2F106X_upgrade2018_realistic_v15_L1v1-v1%2F120000%2F16298DD3-C67F-984F-B396-A8FDA07FE2A2.root&instance=prod/global


		DY 2 Jets
Background_Analysis("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/DYJetsToLL_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/100000/25E5160F-9CB6-104A-8F1C-079DB5FF74CF.root", "Out_DY2.root", 360.4, 1)
		DY MassRange
Background_Analysis("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/110000/090897E2-634D-FC40-A500-F55BE44C83FA.root", "Out_DYM.root", 15810.0, 1)

		W 0 Jet
Background_Analysis("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/WJetsToLNu_0J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/260000/D824EF2C-E82C-8343-A030-A0973EADFE5C.root", "Out_W0.root", 53330.0 , 1)
		W 1 Jet
Background_Analysis("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/WJetsToLNu_1J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/00000/180EAB2C-2584-E143-A6D0-DA375CDF647C.root", "Out_W1.root", 8875.0 , 1)
		W 2 Jet
Background_Analysis("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/WJetsToLNu_2J_TuneCP5_13TeV-amcatnloFXFX-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/00000/925ED0F2-0702-B24F-9CDB-B8A446DBD215.root", "Out_W2.root", 3338.0 , 1)

		SingleTop
Background_Analysis("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/ST_tW_top_5f_NoFullyHadronicDecays_TuneCP5_13TeV_PDFWeights-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/2430000/E6419406-998F-C144-836F-6C03FAE7D6DE.root", "Out_top.root", 32.45 , 1)
		SingleAntiTop
Background_Analysis("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/ST_tW_antitop_5f_NoFullyHadronicDecays_TuneCP5_13TeV_PDFWeights-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/2530000/2D9A9938-5A10-9348-A9C6-004CFB825B74.root", "Out_antitop.root", 32.51 , 1)