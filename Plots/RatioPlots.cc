#include<iostream>
#include<TH1.h>
#include<TFile.h>
#include<array>
#include<THStack.h>
#include<vector>
#include<stdio.h>

#define NFILES 7
#define NMCFILES 6

using namespace std;

void Distributions(){
 string startdir="/afs/cern.ch/user/g/gdamolin/Johan/Temp/";
 TFile *FileData = new TFile((startdir+"Data.root").c_str());

 TFile *FileDB = new TFile((startdir+"Diboson.root").c_str());
 TFile *FileDY = new TFile((startdir+"DY.root").c_str());
 TFile *FileST = new TFile((startdir+"ST.root").c_str());
 TFile *FileTT = new TFile((startdir+"TTbkgd.root").c_str());
 TFile *FileW = new TFile((startdir+"W.root").c_str());
 TFile *FileSig = new TFile((startdir+"Signal.root").c_str());

 array<TFile*, NFILES> Files={FileData,FileDB,FileDY, FileST, FileTT,FileW,FileSig};
 array<TH1D*,7> Mupt;
 array<TH1D*,7> Mueta;
 array<TH1D*,7> Ept;
 array<TH1D*,7> Eeta;
 array<TH1D*,7> LeadingPt;
 array<TH1D*,7> InvMass;

 array<TTree*,NMCFILES> Trees;

 for(int i=0;i<NFILES;i++){
  TTree * tout= static_cast<TTree *>(Files[i]->Get("tout"));
   string dset;
   switch(i){
        case 0: dset="Data";break;
	case 1: dset="DiBoson"; break;
        case 2: dset="DY"; break;
	case 3: dset="ST"; break;
	case 4: dset="TTbkg"; break;
	case 5: dset="Wjeta"; break;
	case 6: dset="TTleptonic"; break;
   }
  char nameA[20]; char nameB[20];char nameC[20]; char nameD[20]; char nameE[20]; char nameF[20];
  sprintf(nameA,"Mupt_%s",dset.c_str()); sprintf(nameB,"Mueta_%s",dset.c_str()); sprintf(nameC,"Ept_%s",dset.c_str());
  sprintf(nameD,"Eeta_%s",dset.c_str()); sprintf(nameE,"2LMass_%s",dset.c_str()); sprintf(nameF,"LeadPt_%s",dset.c_str());
  if(i==0){
   Mupt[i] = static_cast<TH1D *>(Files[i]->Get("h_Muon_pt"));
   Mueta[i] = static_cast<TH1D *>(Files[i]->Get("h_Muon_eta"));
   Ept[i] = static_cast<TH1D *>(Files[i]->Get("h_Electron_pt"));
   Eeta[i] = static_cast<TH1D *>(Files[i]->Get("h_Electron_eta"));
   InvMass[i] = static_cast<TH1D *>(Files[i]->Get("h_Muon_Electron_invariant_mass"));
   LeadingPt[i] = static_cast<TH1D *>(Files[i]->Get("h_leading_lepton_pt"));
   }
  if (i>0){
   Mupt[i] = static_cast<TH1D *>(Files[i]->Get("h_Muon_pt_weighted"));
   Mueta[i] = static_cast<TH1D *>(Files[i]->Get("h_Muon_eta_weighted"));
   Ept[i] = static_cast<TH1D *>(Files[i]->Get("h_Electron_pt_weighted"));
   Eeta[i] = static_cast<TH1D *>(Files[i]->Get("h_Electron_eta_weighted"));
   InvMass[i] = static_cast<TH1D *>(Files[i]->Get("h_Muon_Electron_invariant_mass_weighted"));
   LeadingPt[i] = static_cast<TH1D *>(Files[i]->Get("h_leading_lepton_pt_weighted"));
   }
  Mupt[i]->SetNameTitle(nameA,nameA);
  Mueta[i]->SetNameTitle(nameB,nameB);
  Ept[i]->SetNameTitle(nameC,nameC);
  Eeta[i]->SetNameTitle(nameD,nameD);
  InvMass[i]->SetNameTitle(nameE,nameE);
  LeadingPt[i]->SetNameTitle(nameF,nameF);
  
 }

 cout<<"Blue: signal ev         "<< Mueta[6]->Integral()<<endl;
 cout<<"Green: Dibosons         "<< Mueta[1]->Integral()<<endl;
 cout<<"Black: DY               "<< Mueta[2]->Integral()<<endl;
 cout<<"Red: single top         "<< Mueta[3]->Integral()<<endl;
 cout<<"Orange: ttbkg           "<< Mueta[4]->Integral()<<endl;
 cout<<"Yellow: Wjets           "<< Mueta[5]->Integral()<<endl;

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
        Mueta[0]->SetMarkerSize(2);
        Mueta[0]->SetMarkerStyle(2);
        Mueta[0]->Draw("P SAME");
	TRatioPlot* a1 = new TRatioPlot(hs1,Mueta[0]);
	a1->Draw();
        d1->Update();
        gPad->BuildLegend(0.7,0.7,0.9,0.9);
        gPad->Update();

        TCanvas *d2 = new TCanvas("Eeta");
        THStack *hs2 = new THStack("hs2", "");
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
        hs2->Add(Eeta[1]); hs2->Add(Eeta[2]);
        hs2->Add(Eeta[3]); hs2->Add(Eeta[4]);
        hs2->Add(Eeta[5]); hs2->Add(Eeta[6]);
        hs2->Draw("HIST");
        hs2->GetXaxis()->SetTitle("Eta of Electrons");
        Eeta[0]->SetMarkerColor(kRed);
        Eeta[0]->SetMarkerSize(2);
        Eeta[0]->SetMarkerStyle(2);
        Eeta[0]->Draw("P SAME");
	TRatioPlot* a2 = new TRatioPlot(hs2,Eeta[0]);
	a2->Draw();
        d2->Update();
        gPad->BuildLegend(0.7,0.7,0.9,0.9);
        gPad->Update();
        
        TCanvas *d3 = new TCanvas("Ept");
        THStack *hs3 = new THStack("hs3", "");
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
        hs3->Add(Ept[1]); hs3->Add(Ept[2]);
        hs3->Add(Ept[3]); hs3->Add(Ept[4]);
        hs3->Add(Ept[5]); hs3->Add(Ept[6]);
        hs3->Draw("HIST");
        hs3->GetXaxis()->SetTitle("pt of Electrons");
        Ept[0]->SetMarkerColor(kRed);
        Ept[0]->SetMarkerSize(2);
        Ept[0]->SetMarkerStyle(2);
        Ept[0]->Draw("P SAME");
	TRatioPlot* a3 = new TRatioPlot(hs3,Ept[0]);
	a3->Draw();
        d3->Update();
        gPad->BuildLegend(0.7,0.7,0.9,0.9);
        gPad->Update();

        TCanvas *d4 = new TCanvas("Mupt");
        THStack *hs4 = new THStack("hs4", "");
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
        hs4->Add(Mupt[1]); hs4->Add(Mupt[2]);
        hs4->Add(Mupt[3]); hs4->Add(Mupt[4]);
        hs4->Add(Mupt[5]); hs4->Add(Mupt[6]);
        hs4->Draw("HIST");
        hs4->GetXaxis()->SetTitle("pt of Muons");
        Mupt[0]->SetMarkerColor(kRed);
        Mupt[0]->SetMarkerSize(2);
        Mupt[0]->SetMarkerStyle(2);
        Mupt[0]->Draw("P SAME");
	TRatioPlot* a4 = new TRatioPlot(hs4,Mupt[0]);
	a4->Draw();
        d4->Update();
        gPad->BuildLegend(0.7,0.7,0.9,0.9);
        gPad->Update();

	TCanvas *d5 = new TCanvas("LeadingPt");
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
        hs5->GetXaxis()->SetTitle("pt of Leading Lepton");
        LeadingPt[0]->SetMarkerColor(kRed);
        LeadingPt[0]->SetMarkerSize(2);
        LeadingPt[0]->SetMarkerStyle(2);
        LeadingPt[0]->Draw("P SAME");
	TRatioPlot* a5 = new TRatioPlot(hs5,LeadingPt[0]);
	a5->Draw();
        d5->Update();
        gPad->BuildLegend(0.7,0.7,0.9,0.9);
        gPad->Update();

	TCanvas *d6 = new TCanvas("Invariant Mass");
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
        hs6->GetXaxis()->SetTitle("Inv Mass");
        InvMass[0]->SetMarkerColor(kRed);
        InvMass[0]->SetMarkerSize(2);
        InvMass[0]->SetMarkerStyle(2);
        InvMass[0]->Draw("P SAME");
	TRatioPlot* a6 = new TRatioPlot(hs6,InvMass[0]);
	a6->Draw();
        d6->Update();
        gPad->BuildLegend(0.7,0.7,0.9,0.9);
        gPad->Update();

return;
}
