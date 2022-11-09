#include<iostream>
#include<TH1.h>
#include<TFile.h>


void ormalize(TH1D *a, Double_t norm=1, int nbins = 5)
{
	a->Scale(1./norm);
	//a->Rebin(nbins);
}

void Normalizer(string infile){
 TFile *FileIn = TFile::Open(infile.c_str(),"Read");
 //he complains if the tree is big because it's open in read and he cannot save the buffer branch anywhere

 TTree * runIn= static_cast<TTree *>(FileIn->Get("Run_out"));
 TTree * tin= static_cast<TTree *>(FileIn->Get("tout"));

 Float_t CorrWeight,Weight;
 tin->SetBranchStatus("*","0");
 tin->SetBranchStatus("Weight","1");
 tin->SetBranchAddress("Weight",&Weight);
 
 TBranch *bout = tin->Branch("CorrWeight",&CorrWeight,"CorrWeight/F");


 Double_t singFile=0, TotW=0;
 runIn->SetBranchStatus("*",0);
 runIn->SetBranchStatus("genEventSumw",1);
 runIn->SetBranchAddress("genEventSumw",&singFile);

 for(int j=0; j<runIn->GetEntries();j++){
   	runIn->GetEntry(j);
    	TotW+=singFile;
   	}

 for(int j=0; j<tin->GetEntries();j++){
   	tin->GetEntry(j);
    	CorrWeight=Weight/TotW;
        bout->Fill();
   	}

 TH1D* Mupt = static_cast<TH1D *>(FileIn->Get("h_Muon1_pt"));
 TH1D* Mueta = static_cast<TH1D *>(FileIn->Get("h_Muon1_eta"));
 TH1D* Ept = static_cast<TH1D *>(FileIn->Get("h_Muon2_pt"));
 TH1D* Eeta = static_cast<TH1D *>(FileIn->Get("h_Muon2_eta"));
 TH1D* InvMass = static_cast<TH1D *>(FileIn->Get("Muon_Muon_invariant_mass"));
 TH1D* LeadingPt = static_cast<TH1D *>(FileIn->Get("h_Acopla_emu"));
 TH1D* NJET = static_cast<TH1D *>(FileIn->Get("N_Jets"));

 ormalize(Mupt,TotW); ormalize(Mueta,TotW); ormalize(Ept,TotW); ormalize(Eeta,TotW); ormalize(InvMass,TotW); ormalize(LeadingPt,TotW); ormalize(NJET,TotW);

 std::cout<<"Normalization "<<TotW<<std::endl;

 TFile *FileOut = TFile::Open(("Corr/"+infile).c_str(),"RECREATE");
 TTree *tout=tin->CloneTree();
//FileOut->cd();
 Mupt->Write();
 Mueta->Write();
 Ept->Write();
 Eeta->Write();
 InvMass->Write();
 LeadingPt->Write();
 NJET->Write();
 tout->Write();

 FileOut->Write();
 FileOut->Close();
 return;
 }
 
