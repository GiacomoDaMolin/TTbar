#include<iostream>
#include<TH1.h>
#include<TFile.h>

void Normalize(TH1D *a, Double_t norm=1, int nbins = 5)
{
	a->Scale(1./norm*1000);
	//a->Rebin(nbins);
}

void Normalizer(string infile){
 string startdir="/eos/user/g/gdamolin/TT2bbemu/RESULTS/";
 TFile *FileIn = new TFile((startdir+infile).c_str());

 TTree * runIn= static_cast<TTree *>(FileIn->Get("Run_out"));
 Double_t singFile=0, TotW=0;
 runIn->SetBranchStatus("*",0);
 runIn->SetBranchStatus("genEventSumw",1);
 runIn->SetBranchAddress("genEventSumw",&singFile);
  	for(int j=0; j<runIn->GetEntries();j++){
   	runIn->GetEntry(j);
    	TotW+=singFile;
	//cout<<TotW<<endl;
   	}
 TH1D* Mupt = static_cast<TH1D *>(FileIn->Get("h_Muon_pt_weighted"));
 TH1D* Mueta = static_cast<TH1D *>(FileIn->Get("h_Muon_eta_weighted"));
 TH1D* Ept = static_cast<TH1D *>(FileIn->Get("h_Electron_pt_weighted"));
 TH1D* Eeta = static_cast<TH1D *>(FileIn->Get("h_Electron_eta_weighted"));
 TH1D* InvMass = static_cast<TH1D *>(FileIn->Get("Muon_Electron_invariant_mass_weighted"));
 TH1D* LeadingPt = static_cast<TH1D *>(FileIn->Get("leading_lepton_pt_weighted"));

 Normalize(Mupt,TotW); Normalize(Mueta,TotW); Normalize(Ept,TotW); Normalize(Eeta,TotW); Normalize(InvMass,TotW); Normalize(LeadingPt,TotW);

 TFile *FileOut = new TFile((startdir+"Corr/"+infile).c_str(),"RECREATE");
 Mupt->Write();
 Mueta->Write();
 Ept->Write();
 Eeta->Write();
 InvMass->Write();
 LeadingPt->Write();
 //runIn->Write();

 FileOut->Write();
 FileOut->Close();
 return;
 }
 