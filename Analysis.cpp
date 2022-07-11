#include <iostream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"


using namespace std;

void func(string inputFile, string ofile);

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

    //const auto nEv = tin->GetEntries();
    //for (int i = 0; i < nEv; i++){
    //    
    //}
    
}