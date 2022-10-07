#include<iostream>
#include<TH1.h>
#include<TFile.h>
#include<array>
#include<THStack.h>

void Normalize(TH1D *a, Double_t norm=1, int nbins = 5)
{
	a->Scale(1./norm*1000);
	//a->Rebin(nbins);
}
void Distributions(){
 string startdir="/afs/cern.ch/user/g/gdamolin/public/DR/Corr/";
 TFile *FileData = new TFile((startdir+"Out_Data.root").c_str());

 TFile *FileDB = new TFile((startdir+"Out_DB.root").c_str()); 
 TFile *FileDY = new TFile((startdir+"Out_DY.root").c_str());
 TFile *FileST = new TFile((startdir+"Out_ST.root").c_str());
 TFile *FileTT = new TFile((startdir+"Out_TTbkg.root").c_str());
 TFile *FileW = new TFile((startdir+"Out_W.root").c_str());
 TFile *FileSig = new TFile((startdir+"Out_Signal.root").c_str());
  
 array<TFile*, 7> Files={FileData,FileDB,FileDY, FileST, FileTT,FileW,FileSig};
 array<TH1D*,7> Mupt;
 array<TH1D*,7> Mueta;
 array<TH1D*,7> Ept;
 array<TH1D*,7> Eeta;
 array<TH1D*,7> LeadingPt;
 array<TH1D*,7> InvMass;

 array<TTree*,6> Trees;
//loop in each file, get histograms, normalize if necessary save them in arrays
 //cout<<"Start of the cycle"<<endl;
 for(int i=0;i<7;i++){
  //cout<<"########### File "<< i <<endl;
  if(i==0){
   Mupt[i] = static_cast<TH1D *>(Files[i]->Get("h_Muon_pt"));
   Mueta[i] = static_cast<TH1D *>(Files[i]->Get("h_Muon_eta"));
   Ept[i] = static_cast<TH1D *>(Files[i]->Get("h_Electron_pt"));
   Eeta[i] = static_cast<TH1D *>(Files[i]->Get("h_Electron_eta"));
   InvMass[i] = static_cast<TH1D *>(Files[i]->Get("Muon_Electron_invariant_mass"));
   LeadingPt[i] = static_cast<TH1D *>(Files[i]->Get("leading_lepton_pt"));
   //Normalize(Mupt[i]); Normalize(Mueta[i]); Normalize(Ept[i]); Normalize(Eeta[i]); Normalize(InvMass[i]); Normalize(LeadingPt[i]);
   }
 if(i>0) {
 /* TTree * runout= static_cast<TTree *>(Files[i]->Get("Run_out"));
  Double_t singFile=0, TotW=0;
  runout->SetBranchStatus("*",0);
  runout->SetBranchStatus("genEventSumw",1);
  runout->SetBranchAddress("genEventSumw",&singFile);
  	for(int j=0; j<runout->GetEntries();j++){
   	runout->GetEntry(j);
    	TotW+=singFile;
	//cout<<TotW<<endl;
   	}*/
  Mupt[i] = static_cast<TH1D *>(Files[i]->Get("h_Muon_pt_weighted"));
  Mueta[i] = static_cast<TH1D *>(Files[i]->Get("h_Muon_eta_weighted"));
  Ept[i] = static_cast<TH1D *>(Files[i]->Get("h_Electron_pt_weighted"));
  Eeta[i] = static_cast<TH1D *>(Files[i]->Get("h_Electron_eta_weighted"));
  InvMass[i] = static_cast<TH1D *>(Files[i]->Get("Muon_Electron_invariant_mass_weighted"));
  LeadingPt[i] = static_cast<TH1D *>(Files[i]->Get("leading_lepton_pt_weighted"));
  //Normalize(Mupt[i],TotW); Normalize(Mueta[i],TotW); Normalize(Ept[i],TotW); Normalize(Eeta[i],TotW); Normalize(InvMass[i],TotW); Normalize(LeadingPt[i],TotW);
  }
 } 

//Plot from arrays
 cout<<"Blue: signal ev		"<< Mueta[6]->Integral()<<endl;
 cout<<"Green: Dibosons		"<< Mueta[1]->Integral()<<endl;
 cout<<"Black: DY		"<< Mueta[2]->Integral()<<endl;
 cout<<"Red: single top		"<< Mueta[3]->Integral()<<endl;
 cout<<"Orange: ttbkg		"<< Mueta[4]->Integral()<<endl;
 cout<<"Yellow: Wjets		"<< Mueta[5]->Integral()<<endl;

 cout<<"MCsum "<< Mueta[6]->Integral()+ Mueta[1]->Integral()+ Mueta[2]->Integral()+Mueta[3]->Integral()+ Mueta[4]->Integral()+Mueta[5]->Integral()<<endl;
 cout<<"Data "<< Mueta[0]->Integral()<<endl;
 	TCanvas *d1 = new TCanvas("Mueta");
	THStack *hs1 = new THStack("hs1", "");
	Mueta[6]->SetFillColor(kBlue);
	Mueta[1]->SetFillColor(kGreen);
	Mueta[2]->SetFillColor(kBlack);
	Mueta[3]->SetFillColor(kRed);
	Mueta[4]->SetFillColor(kOrange);
	Mueta[5]->SetFillColor(kYellow);
	Mueta[6]->SetMarkerColor(kBlue);
	Mueta[1]->SetMarkerColor(kGreen);
	Mueta[2]->SetMarkerColor(kBlack);
	Mueta[3]->SetMarkerColor(kRed);
	Mueta[4]->SetMarkerColor(kOrange);
	Mueta[5]->SetMarkerColor(kYellow);
	hs1->Add(Mueta[1]); hs1->Add(Mueta[2]);
	hs1->Add(Mueta[3]); hs1->Add(Mueta[4]);
	hs1->Add(Mueta[5]); hs1->Add(Mueta[6]);
	hs1->Draw("HIST");
	hs1->GetXaxis()->SetTitle("Eta of muons");
	Mueta[0]->SetMarkerColor(kRed);
	/*Mueta[0]->SetLineWidth(3);
	Mueta[0]->SetLineColor(kRed);*/
	Mueta[0]->SetMarkerSize(2);
	Mueta[0]->SetMarkerStyle(2);
	Mueta[0]->Draw("P SAME");
	gPad->BuildLegend(0.7,0.7,0.9,0.9);
	gPad->Update();
	d1->Update(); 

	TCanvas *d2 = new TCanvas("mupt");
	THStack *hs2 = new THStack("hs2", "");
	Mupt[6]->SetFillColor(kBlue);
	Mupt[1]->SetFillColor(kGreen);
	Mupt[2]->SetFillColor(kBlack);
	Mupt[3]->SetFillColor(kRed);
	Mupt[4]->SetFillColor(kOrange);
	Mupt[5]->SetFillColor(kYellow);
	Mupt[6]->SetMarkerColor(kBlue);
	Mupt[1]->SetMarkerColor(kGreen);
	Mupt[2]->SetMarkerColor(kBlack);
	Mupt[3]->SetMarkerColor(kRed);
	Mupt[4]->SetMarkerColor(kOrange);
	Mupt[5]->SetMarkerColor(kYellow);
	hs2->Add(Mupt[1]); hs2->Add(Mupt[2]);
	hs2->Add(Mupt[3]); hs2->Add(Mupt[4]);
	hs2->Add(Mupt[5]); hs2->Add(Mupt[6]);
	hs2->Draw("HIST");
	hs2->GetXaxis()->SetTitle("PT of muons");
	Mupt[0]->SetMarkerColor(kRed);
	/*Mupt[0]->SetLineWidth(3);
	Mupt[0]->SetLineColor(kRed);*/
	Mupt[0]->SetMarkerSize(2);
	Mupt[0]->SetMarkerStyle(2);
	Mupt[0]->Draw("P SAME");
	gPad->BuildLegend(0.7,0.7,0.9,0.9);
	gPad->Update();
	d2->Update(); 

	TCanvas *d3 = new TCanvas("EPT");
	THStack *hs3 = new THStack("hs3", "");
	Eeta[6]->SetFillColor(kBlue);
	Eeta[1]->SetFillColor(kGreen);
	Eeta[2]->SetFillColor(kBlack);
	Eeta[3]->SetFillColor(kRed);
	Eeta[4]->SetFillColor(kOrange);
	Eeta[5]->SetFillColor(kYellow);
	Eeta[6]->SetMarkerColor(kBlue);
	Eeta[1]->SetMarkerColor(kGreen);
	Eeta[2]->SetMarkerColor(kBlack);
	Eeta[3]->SetMarkerColor(kRed);
	Eeta[4]->SetMarkerColor(kOrange);
	Eeta[5]->SetMarkerColor(kYellow);
	hs3->Add(Eeta[1]); hs3->Add(Eeta[2]);
	hs3->Add(Eeta[3]); hs3->Add(Eeta[4]);
	hs3->Add(Eeta[5]); hs3->Add(Eeta[6]);
	hs3->Draw("HIST");
	hs3->GetXaxis()->SetTitle("Eta of E");
	Eeta[0]->SetMarkerColor(kRed);
	/*Eeta[0]->SetLineWidth(3);
	Eeta[0]->SetLineColor(kRed);*/
	Eeta[0]->SetMarkerSize(2);
	Eeta[0]->SetMarkerStyle(2);
	Eeta[0]->Draw("P SAME");
	gPad->BuildLegend(0.7,0.7,0.9,0.9);
	gPad->Update();
	d3->Update(); 

	TCanvas *d4 = new TCanvas("Emu");
	THStack *hs4 = new THStack("hs4", "");
	Ept[6]->SetFillColor(kBlue);
	Ept[1]->SetFillColor(kGreen);
	Ept[2]->SetFillColor(kBlack);
	Ept[3]->SetFillColor(kRed);
	Ept[4]->SetFillColor(kOrange);
	Ept[5]->SetFillColor(kYellow);
	Ept[6]->SetMarkerColor(kBlue);
	Ept[1]->SetMarkerColor(kGreen);
	Ept[2]->SetMarkerColor(kBlack);
	Ept[3]->SetMarkerColor(kRed);
	Ept[4]->SetMarkerColor(kOrange);
	Ept[5]->SetMarkerColor(kYellow);
	hs4->Add(Ept[1]); hs4->Add(Ept[2]);
	hs4->Add(Ept[3]); hs4->Add(Ept[4]);
	hs4->Add(Ept[5]); hs4->Add(Ept[6]);
	hs4->Draw("HIST");
	hs4->GetXaxis()->SetTitle("PT of E");
	Ept[0]->SetMarkerColor(kRed);
	/*Ept[0]->SetLineWidth(3);
	Ept[0]->SetLineColor(kRed);*/
	Ept[0]->SetMarkerSize(2);
	Ept[0]->SetMarkerStyle(2);
	Ept[0]->Draw("P SAME");
	gPad->BuildLegend(0.7,0.7,0.9,0.9);
	gPad->Update();
	d4->Update();

	TCanvas *d5 = new TCanvas("lead");
	THStack *hs5 = new THStack("hs5", "");
	LeadingPt[6]->SetFillColor(kBlue);
	LeadingPt[1]->SetFillColor(kGreen);
	LeadingPt[2]->SetFillColor(kBlack);
	LeadingPt[3]->SetFillColor(kRed);
	LeadingPt[4]->SetFillColor(kOrange);
	LeadingPt[5]->SetFillColor(kYellow);
	LeadingPt[6]->SetMarkerColor(kBlue);
	LeadingPt[1]->SetMarkerColor(kGreen);
	LeadingPt[2]->SetMarkerColor(kBlack);
	LeadingPt[3]->SetMarkerColor(kRed);
	LeadingPt[4]->SetMarkerColor(kOrange);
	LeadingPt[5]->SetMarkerColor(kYellow);
	hs5->Add(LeadingPt[1]); hs5->Add(LeadingPt[2]);
	hs5->Add(LeadingPt[3]); hs5->Add(LeadingPt[4]);
	hs5->Add(LeadingPt[5]); hs5->Add(LeadingPt[6]);
	hs5->Draw("HIST");
	hs5->GetXaxis()->SetTitle("PT of Leading lepton [GeV]");
	LeadingPt[0]->SetMarkerColor(kRed);
	/*LeadingPt[0]->SetLineWidth(3);
	LeadingPt[0]->SetLineColor(kRed);*/
	LeadingPt[0]->SetMarkerSize(2);
	LeadingPt[0]->SetMarkerStyle(2);
	LeadingPt[0]->Draw("P SAME");
	gPad->BuildLegend(0.7,0.7,0.9,0.9);
	gPad->Update();
	d5->Update(); 

	TCanvas *d6 = new TCanvas("maas");
	THStack *hs6 = new THStack("hs6", "");
	InvMass[6]->SetFillColor(kBlue);
	InvMass[1]->SetFillColor(kGreen);
	InvMass[2]->SetFillColor(kBlack);
	InvMass[3]->SetFillColor(kRed);
	InvMass[4]->SetFillColor(kOrange);
	InvMass[5]->SetFillColor(kYellow);
	InvMass[6]->SetMarkerColor(kBlue);
	InvMass[1]->SetMarkerColor(kGreen);
	InvMass[2]->SetMarkerColor(kBlack);
	InvMass[3]->SetMarkerColor(kRed);
	InvMass[4]->SetMarkerColor(kOrange);
	InvMass[5]->SetMarkerColor(kYellow);
	hs6->Add(InvMass[1]); hs6->Add(InvMass[2]);
	hs6->Add(InvMass[3]); hs6->Add(InvMass[4]);
	hs6->Add(InvMass[5]); hs6->Add(InvMass[6]);
	hs6->Draw("HIST");
	hs6->GetXaxis()->SetTitle("Invariant mass of leptonic system [GeV]");
	InvMass[0]->SetMarkerColor(kRed);
	/*InvMass[0]->SetLineWidth(3);
	InvMass[0]->SetLineColor(kRed);*/
	InvMass[0]->SetMarkerSize(2);
	InvMass[0]->SetMarkerStyle(2);
	InvMass[0]->Draw("P SAME");
	gPad->BuildLegend(0.7,0.7,0.9,0.9);
	gPad->Update();
	d6->Update();

}
