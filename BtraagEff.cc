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
//TFile *fin= TFile::Open("root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18NanoAODv2/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/NANOAODSIM/106X_upgrade2018_realistic_v15_L1v1-v1/270000/A62BD1E4-52E1-F64F-9F99-894D7A4B19D1.root");
 TTree *tin = static_cast<TTree *>(fin->Get("Events"));
 tin->SetBranchStatus("*", 0);
 Float_t Jet_eta[MAX_ARRAY_SIZE],Jet_btagDeepFlavB[MAX_ARRAY_SIZE],Jet_pt[MAX_ARRAY_SIZE];
 UInt_t nJet;
 Int_t Jet_jetId[MAX_ARRAY_SIZE], Jet_puId[MAX_ARRAY_SIZE];
 Int_t Jet_hadronFlavour[MAX_ARRAY_SIZE];
 std::cout<<"Set branches"<<std::endl;
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
 
 double ptbins[9]={25,40,60,80,100,140,180,260,1000};

 double etabins[5]={0, 0.6, 1.2, 1.8, 2.4};

 std::cout<<"Create histos"<<std::endl; 
 //make unequal size bins: you have to have a bin for every pT. so a very large bin for very large pT is needed
 TH2D * h2_BTaggingEff_Denom_b= new TH2D("b_jets","b_jets vs pt and eta; p_T [GeV]; Pseudorapidity",8,ptbins,4,etabins);
 TH2D * h2_BTaggingEff_Denom_c= new TH2D("c_jets","c_jets vs pt and eta; p_T [GeV]; Pseudorapidity",8,ptbins,4,etabins);
 TH2D * h2_BTaggingEff_Denom_udsg= new TH2D("l_jets","l_jets vs pt and eta; p_T [GeV]; Pseudorapidity",8,ptbins,4,etabins);
 TH2D * h2_BTaggingEff_Num_b= new TH2D("b_jets_tagged","Tagged b_jets vs pt and eta; p_T [GeV]; Pseudorapidity",8,ptbins,4,etabins);
 TH2D * h2_BTaggingEff_Num_c= new TH2D("c_jets_tagged","Tagged c_jets vs pt and eta; p_T [GeV]; Pseudorapidity",8,ptbins,4,etabins);
 TH2D * h2_BTaggingEff_Num_udsg= new TH2D("l_jets_tagged","Tagged light_jets vs pt and eta; p_T [GeV]; Pseudorapidity",8,ptbins,4,etabins);

 h2_BTaggingEff_Denom_b->Sumw2();
 h2_BTaggingEff_Denom_c->Sumw2();
 h2_BTaggingEff_Denom_udsg->Sumw2();
 h2_BTaggingEff_Num_b->Sumw2();
 h2_BTaggingEff_Num_c->Sumw2();
 h2_BTaggingEff_Num_udsg->Sumw2();
 

 //#pragma omp parallel for
 for(int i=0; i<tin->GetEntries(); i++){
  tin->GetEntry(i);
  if (i % 100000 == 0)
            std::cout << "Processing entry " << i << " of " << tin->GetEntries() << std::endl;
  for (int j=0; j<nJet;j++){
    if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6) && (Jet_pt[j]>50 || Jet_puId[j]==7)){
    int flav=Jet_hadronFlavour[j];
    /*if((abs(abs(Jet_eta[j])-1.5)<0.3)) {std::cout<<Jet_eta[j]<<std::endl;
	std::cout<<abs(abs(Jet_eta[j])-1.5)<<std::endl;}*/
    if(flav==5) h2_BTaggingEff_Denom_b->Fill(Jet_pt[j],abs(Jet_eta[j]));
    if(flav==4) h2_BTaggingEff_Denom_c->Fill(Jet_pt[j],abs(Jet_eta[j]));
    if(flav==0) h2_BTaggingEff_Denom_udsg->Fill(Jet_pt[j],abs(Jet_eta[j]));
    if(Jet_btagDeepFlavB[j]>0.2783){
    if(flav==5) h2_BTaggingEff_Num_b->Fill(Jet_pt[j],abs(Jet_eta[j]));
    if(flav==4) h2_BTaggingEff_Num_c->Fill(Jet_pt[j],abs(Jet_eta[j]));
    if(flav==0) h2_BTaggingEff_Num_udsg->Fill(Jet_pt[j],abs(Jet_eta[j]));
    }
   } 
  }
 }
 auto const pos = infile.find_last_of('/');
 std::string inname=infile.substr(pos + 1);
 std::string outfile="/afs/cern.ch/user/g/gdamolin/public/tempB/Out_"+inname; //DO NOT DO HERE THE RATIO, YOU CANNOT HADD FILES IF YOU DO
 std::cout<<outfile<<std::endl;
 TFile *a=TFile::Open(outfile.c_str(),"recreate");
 h2_BTaggingEff_Denom_b->Write();
 h2_BTaggingEff_Denom_c->Write();
 h2_BTaggingEff_Denom_udsg->Write();
 h2_BTaggingEff_Num_b->Write();
 h2_BTaggingEff_Num_c->Write();
 h2_BTaggingEff_Num_udsg->Write();
 
 a->Close();
}


int main(int argc, char **argv)
{

    std::string inputFile = argv[1];
    std::cout << inputFile.c_str() <<std::endl;
    EffComputer(inputFile);

    return 0;
}
