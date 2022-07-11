#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"


using namespace std;

#define MAX_ARRAY_SIZE 1024

int main(int argc, char** argv){
    string inputFile = argv[1];
    string outputFile = argv[2];
    func(inputFile, outputFile);
}
void func(string inputFile, string ofile){

    TFile *fin = TFile::Open(inputFile.c_str());
    TTree *tin = static_cast<TTree*>(fin->Get("Events"));

    TH1D* myhist = new TH1D("PTmuon","#\\mu p_{\\mathrm{T}}", 150, 0,150);
    tin->Draw("Muon_pt >> PTmuon", "(HLT_IsoMu27) ", "GOFF");

    // myhist.Draw() will automatically draw to the current 'workspace' which is mycanvas
    TCanvas* mycanvas = new TCanvas("MuonCanvas", "MuonCanvas", 800, 600);
    myhist->Draw();
    mycanvas->Draw();
    mycanvas->Update();
    int var;

    Float_t Muon_pt[MAX_ARRAY_SIZE], Electron_pt[MAX_ARRAY_SIZE];
    Float_t Muon_eta[MAX_ARRAY_SIZE], Electron_eta[MAX_ARRAY_SIZE];
    Float_t Muon_phi[MAX_ARRAY_SIZE], Electron_phi[MAX_ARRAY_SIZE];
    Float_t Muon_mass[MAX_ARRAY_SIZE], Electron_mass[MAX_ARRAY_SIZE];
    UInt_t nMuon, nElectron;
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



    const auto nEv = tin->GetEntries();
    for (size_t i = 0; i < nEv; i++){
        tin->GetEntry(i)
        
    }
    
}