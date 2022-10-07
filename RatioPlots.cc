#include<iostream>
#include<TH1.h>
#include<TFile.h>
#include<array>
#include<THStack.h>

#define NFILES 7
#define NMCFILES 6

void Distributions(){
 string startdir="/afs/cern.ch/user/g/gdamolin/public/DR/Corr/";
 TFile *FileData = new TFile((startdir+"Out_Data.root").c_str());

 TFile *FileDB = new TFile((startdir+"Out_DB.root").c_str());
 TFile *FileDY = new TFile((startdir+"Out_DY.root").c_str());
 TFile *FileST = new TFile((startdir+"Out_ST.root").c_str());
 TFile *FileTT = new TFile((startdir+"Out_TTbkg.root").c_str());
 TFile *FileW = new TFile((startdir+"Out_W.root").c_str());
 TFile *FileSig = new TFile((startdir+"Out_Signal.root").c_str());

 array<TFile*, NFILES> Files={FileData,FileDB,FileDY, FileST, FileTT,FileW,FileSig};
 array<TH1D*,NFILES> Mupt;
 array<TH1D*,NFILES> Mueta;
 array<TH1D*,NFILES> Ept;
 array<TH1D* NFILES> Eeta;
 array<TH1D*,NFILES> LeadingPt;
 array<TH1D*,NFILES> InvMass;

 array<TTree*,NMCFILES> Trees;

 for(int i=0;i<NFILES;i++){
  //cout<<"########### File "<< i <<endl;
  TTree * tout= static_cast<TTree *>(Files[i]->Get("tout"));
  if(i==0){
   tout->Draw("Muon_pt>>Mupt[i]","",“goff”);
   tout->Draw("Muon_eta>>Mueta[i]","",“goff”);
   tout->Draw("Electron_pt>>Ept[i]","",“goff”);
   tout->Draw("Electron_eta>>Eeta[i]","",“goff”);
   tout->Draw("mu_e_inv_mass>>InvMass[i]","",“goff”);
   tout->Draw("leading_lepton_pt>>LeadingPt[i]","",“goff”);
   }
 if(i>0) {
  tout->Draw("Muon_pt>>Mupt[i]","norm_W>>Mupt[i]",“goff”);
  tout->Draw("Muon_eta>>Mueta[i]","norm_W>>Mueta[i]",“goff”);
  tout->Draw("Electron_pt>>Ept[i]","norm_W>>Ept[i]",“goff”);
  tout->Draw("Electron_eta>>Eeta[i]","norm_W>>Eeta[i]",“goff”);
  tout->Draw("mu_e_inv_mass>>InvMass[i]","norm_W>>InvMass[i]",“goff”);
  tout->Draw("leading_lepton_pt>>LeadingPt[i]","norm_W>>LeadingPt[i]",“goff”);
  }
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
        gPad->BuildLegend(0.7,0.7,0.9,0.9);
        gPad->Update();
        d1->Update();



return;
}




 return;
}
