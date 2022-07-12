#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"


using namespace std;

#define MAX_ARRAY_SIZE 64


void func(string inputFile, string ofile);
int main(int argc, char** argv){

  /*array<TTree*,N> TinArr
    for (int i=0; i++;i<argc){
	TFile a=TFile(argv[1]);
	TTree *ti = new TTree("Tree"+argv[1]");
	TinArr[i]=ti;
	}
    */
    string inputFile = argv[1];
    string outputFile = argv[2];
    func(inputFile, outputFile);
}

void func(string inputFile, string ofile){

    TFile *fin = TFile::Open(inputFile.c_str());
    TTree *tin = static_cast<TTree*>(fin->Get("Events"));

    //TH1D* myhist = new TH1D("PTmuon","#\\mu p_{\\mathrm{T}}", 150, 0,150);
    //tin->Draw("Muon_pt >> PTmuon", "(HLT_IsoMu27) ", "GOFF");

    // myhist.Draw() will automatically draw to the current 'workspace' which is mycanvas
    //TCanvas* mycanvas = new TCanvas("MuonCanvas", "MuonCanvas", 800, 600);
    //myhist->Draw();
    //mycanvas->Draw();
    //mycanvas->Update();
    //int var;

    Float_t Muon_pt[MAX_ARRAY_SIZE], Electron_pt[MAX_ARRAY_SIZE], Jet_pt[MAX_ARRAY_SIZE];
    Float_t Muon_eta[MAX_ARRAY_SIZE], Electron_eta[MAX_ARRAY_SIZE], Jet_eta[MAX_ARRAY_SIZE];
    Float_t Muon_phi[MAX_ARRAY_SIZE], Electron_phi[MAX_ARRAY_SIZE], Jet_phi[MAX_ARRAY_SIZE];
    Float_t Muon_mass[MAX_ARRAY_SIZE], Electron_mass[MAX_ARRAY_SIZE], Jet_mass[MAX_ARRAY_SIZE];
    UInt_t nMuon, nElectron;
    // Set all branches to 0
    tin->SetBranchStatus("*", 0);
    // get the pt
    tin->SetBranchStatus("Electron_pt", 1);
    tin->SetBranchAddress("Electron_pt", &Electron_pt);
    tin->SetBranchStatus("Muon_pt", 1);
    tin->SetBranchAddress("Muon_pt", &Muon_pt);
    // get the number of muons, electrons
    tin->SetBranchStatus("nElectron", 1);
    tin->SetBranchAddress("nElectron", &nElectron);
    tin->SetBranchStatus("nMuon", 1);
    tin->SetBranchAddress("nMuon", &nMuon);
    // get the eta
    tin->SetBranchStatus("Electron_eta", 1);
    tin->SetBranchAddress("Electron_eta", &Electron_eta);
    tin->SetBranchStatus("Muon_eta", 1);
    tin->SetBranchAddress("Muon_eta", &Muon_eta);
    // get the phi
    tin->SetBranchStatus("Electron_phi", 1);
    tin->SetBranchAddress("Electron_phi", &Electron_phi);
    tin->SetBranchStatus("Muon_phi", 1);
    tin->SetBranchAddress("Muon_phi", &Muon_phi);
    // get the mass
    tin->SetBranchStatus("Electron_mass", 1);
    tin->SetBranchAddress("Electron_mass", &Electron_mass);
    tin->SetBranchStatus("Muon_mass", 1);
    tin->SetBranchAddress("Muon_mass", &Muon_mass);
    // get jet quantities
    tin->SetBranchStatus("Jet_pt", 1);
    tin->SetBranchAddress("Jet_pt", &Jet_pt);
    tin->SetBranchStatus("Jet_eta", 1);
    tin->SetBranchAddress("Jet_eta", &Jet_eta);
    tin->SetBranchStatus("Jet_phi", 1);
    tin->SetBranchAddress("Jet_phi", &Jet_phi);
    tin->SetBranchStatus("Jet_mass", 1);
    tin->SetBranchAddress("Jet_mass", &Jet_mass);

    // collect the trigger information
    Bool_t HLT_IsoMu27, HLT_Ele35_WPTight_Gsf; 
    tin->SetBranchStatus("HLT_IsoMu27", 1);
    tin->SetBranchStatus("HLT_Ele35_WPTight_Gsf", 1);
    tin->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27);
    tin->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf);


    // Define the histogram objects
    TH1D* h_Muon_pt = new TH1D("Muon_pt", "Muon_pt", 100, 0, 200);
    TH1D* h_Muon_eta = new TH1D("Muon_eta", "Muon_eta", 100, -5, 5);
    TH1D* h_Electron_pt = new TH1D("Electron_pt", "Electron_pt", 100, 0, 200);
    TH1D* h_Electron_eta = new TH1D("Electron_eta", "Electron_eta", 100, -5, 5);

    TH1D* h_Muon_pt_trigger = new TH1D("Muon_pt_trigger", "Muon_pt_trigger", 100, 0, 200);
    TH1D* h_Muon_eta_trigger = new TH1D("Muon_eta_trigger", "Muon_eta_trigger", 100, -5, 5);
    TH1D* h_Electron_pt_trigger = new TH1D("Electron_pt_trigger", "Electron_pt_trigger", 100, 0, 200);
    TH1D* h_Electron_eta_trigger = new TH1D("Electron_eta_trigger", "Electron_eta_trigger", 100, -5, 5);


    //tin->Draw("Muon_pt >> h_Muon_pt");
    /*tin->Draw("Muon_eta >> h_Muon_eta","", "GOFF");
    tin->Draw("Electron_pt >> h_Electron_pt","", "GOFF");
    tin->Draw("Electron_eta >> h_Electron_eta","", "GOFF");*/
    //tin->Draw("Muon_pt >> h_Muon_pt_trigger","(HLT_IsoMu27) || (HLT_Ele35_WPTight_Gsf)", "GOFF");
    //tin->Draw("Muon_eta >> h_Muon_eta_trigger", "(HLT_IsoMu27) || (HLT_Ele35_WPTight_Gsf)","GOFF");
    //tin->Draw("Electron_pt >> h_Electron_pt_trigger", "(HLT_IsoMu27) || (HLT_Ele35_WPTight_Gsf)","GOFF");
    //tin->Draw("Electron_eta >> h_Electron_eta_trigger", "(HLT_IsoMu27) || (HLT_Ele35_WPTight_Gsf)","GOFF");
    const auto nEv = tin->GetEntries();
    for (size_t i = 0; i < nEv; i++){
        tin->GetEntry(i);
        for (int j = 0; j<nMuon; j++){
            h_Muon_pt->Fill(Muon_pt[j]);
            h_Muon_eta->Fill(Muon_eta[j]);

            if (HLT_IsoMu27 || HLT_Ele35_WPTight_Gsf){
                h_Muon_pt_trigger->Fill(Muon_pt[j]);
                h_Muon_eta_trigger->Fill(Muon_eta[j]);
            }
        }
        
        for (int j = 0; j<nElectron; j++){
            h_Electron_pt->Fill(Electron_pt[j]);
            h_Electron_eta->Fill(Electron_eta[j]);

            if (HLT_IsoMu27 || HLT_Ele35_WPTight_Gsf){
                h_Electron_pt_trigger->Fill(Electron_pt[j]);
                h_Electron_eta_trigger->Fill(Electron_eta[j]);
            }
        }
    }
    
    
    
    //save the histograms in a new File
    TFile *fout =new TFile(ofile.c_str(),"RECREATE");
    // Write the histograms to the file
    h_Muon_pt->Write(); 
    h_Muon_eta->Write(); 
    h_Electron_eta->Write(); 
    h_Electron_pt->Write();
    h_Muon_pt_trigger->Write();
    h_Muon_eta_trigger->Write();
    h_Electron_pt_trigger->Write();
    h_Electron_eta_trigger->Write();

    fout->Write();
    fout->Close();
    
}
