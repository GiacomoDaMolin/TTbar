#include<TH2F.h>
#include<TFile.h>
#include<TTree.h>
#include<string>
#include<iostream>

#define MAX_ARRAY_SIZE 128


void EffComputer(std::string infile){

 
 TFile *fin= new TFile(infile.c_str());
//TFile *fin= TFile::Open("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/270000/A62BD1E4-52E1-F64F-9F99-894D7A4B19D1.root");
 TTree *tin = static_cast<TTree *>(fin->Get("Events"));
 tin->SetBranchStatus("*", 0);
 Float_t Jet_eta[MAX_ARRAY_SIZE],Jet_btagDeepFlavB[MAX_ARRAY_SIZE],Jet_pt[MAX_ARRAY_SIZE];
 UInt_t nJet;
 Int_t Jet_jetId[MAX_ARRAY_SIZE], Jet_puId[MAX_ARRAY_SIZE];
 Int_t Jet_hadronFlavour[MAX_ARRAY_SIZE];
 
 tin->SetBranchStatus("Jet_pt", 1);
 tin->SetBranchAddress("Jet_pt", &Jet_pt);
 tin->SetBranchStatus("Jet_eta", 1);
 tin->SetBranchAddress("Jet_eta", &Jet_eta);
 tin->SetBranchStatus("Jet_btagDeepFlavB", 1);
 tin->SetBranchAddress("Jet_btagDeepFlavB",&Jet_btagDeepFlavB);
 tin->SetBranchStatus("nJet", 1);
 tin->SetBranchAddress("nJet",&nJet);
 tin->SetBranchStatus("Jet_jetId", 1);
 tin->SetBranchStatus("Jet_puId", 1);
 tin->SetBranchAddress("Jet_jetId", &Jet_jetId);
 tin->SetBranchAddress("Jet_puId", &Jet_puId);
 tin->SetBranchStatus("Jet_hadronFlavour", 1);
 tin->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour);
 
 //make unequal size bins: you have to have a bin for every pT. so a very large bin for very large pT is needed
 TH2D  *h2_BTaggingEff_Denom_b= new TH2D("b_jets","b_jets vs pt and eta",//rest of the constructor);
 TH2D  *h2_BTaggingEff_Denom_c;
 TH2D  *h2_BTaggingEff_Denom_udsg;
 TH2D  *h2_BTaggingEff_Num_b;
 TH2D  *h2_BTaggingEff_Num_c;
 TH2D  *h2_BTaggingEff_Num_udsg;

 

 //#pragma omp parallel for
 for(int i=0; i<tin->GetEntries(); i++){
  tin->GetEntry(i);
  for (UInt_t j=0; j<nJet;j++){
    if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6) && (Jet_pt[j]>50 || Jet_puId[j]==7)){
    int flav=Jet_hadronFlavour[j];
    if(flav==5) h2_BTaggingEff_Denom_b->Fill(Jet_pt[j],Jet_eta[j]);
    if(flav==4) h2_BTaggingEff_Denom_c->Fill(Jet_pt[j],Jet_eta[j]);
    if(flav==0) h2_BTaggingEff_Denom_udsg->Fill(Jet_pt[j],Jet_eta[j]);
    if(Jet_btagDeepFlavB[j]>0.2783){
    if(flav==5) h2_BTaggingEff_Num_b->Fill(Jet_pt[j],Jet_eta[j]);
    if(flav==4) h2_BTaggingEff_Num_c->Fill(Jet_pt[j],Jet_eta[j]);
    if(flav==0) h2_BTaggingEff_Num_udsg->Fill(Jet_pt[j],Jet_eta[j]);
    }
   } 
  }
 }

 std::string outfile="/eos/user/g/gdamolin/EFF/b-tag/Out_"+infile; //DO NOT DO HERE THE RATIO, YOU CANNOT HADD FILES IF YOU DO

 TFile *a=new TFile(outfile.c_str(),"Write");
 h2_BTaggingEff_Denom_b->Write();
 h2_BTaggingEff_Denom_c->Write();
 h2_BTaggingEff_Denom_udsg->Write();
 h2_BTaggingEff_Num_b->Write();
 h2_BTaggingEff_Num_c->Write();
 h2_BTaggingEff_Num_udsg->Write();
 
 a->Write(); a->Close();
}


int main(int argc, char **argv)
{

    std::string inputFile = argv[1];
    
    EffComputer(inputFile);

    return 0;
}
