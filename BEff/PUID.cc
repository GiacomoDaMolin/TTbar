#include<TH2F.h>
#include<TFile.h>
#include<TTree.h>
#include<string>
#include<iostream>

#define MAX_ARRAY_SIZE 128

using std::abs;

void EffComputer(std::string infile){
 std::cout<<"Starting the function"<<std::endl;
 TFile *fin= TFile::Open(("root://cms-xrd-global.cern.ch/"+infile).c_str());
 TTree *tin = static_cast<TTree *>(fin->Get("Events"));
 tin->SetBranchStatus("*", 0);
 Float_t Jet_eta[MAX_ARRAY_SIZE],Jet_pt[MAX_ARRAY_SIZE];
 UInt_t nJet;
 Int_t Jet_jetId[MAX_ARRAY_SIZE], Jet_puId[MAX_ARRAY_SIZE];
 Int_t Jet_genJetIdx[MAX_ARRAY_SIZE];
 std::cout<<"Set branches"<<std::endl;
 tin->SetBranchStatus("Jet_pt", 1);
 tin->SetBranchAddress("Jet_pt", &Jet_pt);
 tin->SetBranchStatus("Jet_eta", 1);
 tin->SetBranchAddress("Jet_eta", &Jet_eta);
 tin->SetBranchStatus("nJet", 1);
 tin->SetBranchAddress("nJet",&nJet);
 tin->SetBranchStatus("Jet_jetId", 1);
 tin->SetBranchStatus("Jet_puId", 1);
 tin->SetBranchAddress("Jet_jetId", &Jet_jetId);
 tin->SetBranchAddress("Jet_puId", &Jet_puId);
 tin->SetBranchStatus("Jet_genJetIdx", 1);
 tin->SetBranchAddress("Jet_genJetIdx", &Jet_genJetIdx);
 
 double ptbins[6]={25,30,35,40,45,50};

 double etabins[9]={-2.4,-1.8,-1.2,-0.6, 0, 0.6, 1.2, 1.8, 2.4};

 std::cout<<"Create histos"<<std::endl; 
 //make unequal size bins: you have to have a bin for every pT. so a very large bin for very large pT is needed
 TH2D * h2_PV_Denom= new TH2D("PV_jets","PV_jets vs pt and eta; p_T [GeV]; Pseudorapidity",5,ptbins,8,etabins);
 TH2D * h2_PU_Denom= new TH2D("PU_jets","PU_jets vs pt and eta; p_T [GeV]; Pseudorapidity",5,ptbins,8,etabins);
 TH2D * h2_PV_Num= new TH2D("PV_jets_tagged","PV_jets_tagged vs pt and eta; p_T [GeV]; Pseudorapidity",5,ptbins,8,etabins);
 TH2D * h2_PU_Num= new TH2D("PU_jets_tagged","PU_jets_tagged vs pt and eta; p_T [GeV]; Pseudorapidity",5,ptbins,8,etabins);

 h2_PV_Denom->Sumw2();
 h2_PU_Denom->Sumw2();
 h2_PV_Num->Sumw2();
 h2_PU_Num->Sumw2();


 #pragma omp parallel for
 for(int i=0; i<tin->GetEntries(); i++){
  tin->GetEntry(i);
  if (i % 100000 == 0)
            std::cout << "Processing entry " << i << " of " << tin->GetEntries() << std::endl;
  for (int j=0; j<nJet;j++){
    if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6) && (Jet_pt[j]<50 )){
    if(Jet_genJetIdx[j]<0) h2_PU_Denom->Fill(Jet_pt[j],Jet_eta[j]);
    if(Jet_genJetIdx[j]>=0) h2_PV_Denom->Fill(Jet_pt[j],Jet_eta[j]);

    if(Jet_puId[j]==7){
    if(Jet_genJetIdx[j]<0) h2_PU_Num->Fill(Jet_pt[j],Jet_eta[j]);
    if(Jet_genJetIdx[j]>=0) h2_PV_Num->Fill(Jet_pt[j],Jet_eta[j]);

    }
   } 
  }
 }
 auto const pos = infile.find_last_of('/');
 std::string inname=infile.substr(pos + 1);
 std::string outfile="/afs/cern.ch/user/g/gdamolin/public/tempB/Out_"+inname; //DO NOT DO HERE THE RATIO, YOU CANNOT HADD FILES IF YOU DO
 std::cout<<outfile<<std::endl;
 TFile *a=TFile::Open(outfile.c_str(),"recreate");
 h2_PV_Denom->Write();
 h2_PU_Denom->Write();
 h2_PV_Num->Write();
 h2_PU_Num->Write();
 
 a->Close();
}


int main(int argc, char **argv)
{

    std::string inputFile = argv[1];
    std::cout << inputFile.c_str() <<std::endl;
    EffComputer(inputFile);

    return 0;
}
