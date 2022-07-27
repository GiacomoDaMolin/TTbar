#include <iostream>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"

void Rebin(string inputFile, string OutputFile, int nbins = 10)
{

    TFile *ifile = TFile::Open(inputFile.c_str());
    TFile *ofile = TFile::Open(OutputFile.c_str(), "RECREATE");
    TH1D *Muon_Electron_invariant_mass = (TH1D*)ifile->Get("Muon_Electron_invariant_mass");
    TH1D *h_Muon_pt = (TH1D*)ifile->Get("h_Muon_pt");
    TH1D *h_Muon_eta = (TH1D*)ifile->Get("h_Muon_eta");
    TH1D *h_Electron_pt = (TH1D*)ifile->Get("h_Electron_pt");
    TH1D *h_Electron_eta = (TH1D*)ifile->Get("h_Electron_eta");
    TH1D *leading_lepton_pt = (TH1D*)ifile->Get("leading_lepton_pt");
    Muon_Electron_invariant_mass->Rebin(nbins);
    Muon_Electron_invariant_mass->Scale(1/1000.);
    Muon_Electron_invariant_mass->Write();
    h_Muon_pt->Rebin(nbins);
    h_Muon_pt->Scale(1/1000.);
    h_Muon_pt->Write();
    h_Muon_eta->Rebin(nbins);
    h_Muon_eta->Scale(1/1000.);
    h_Muon_eta->Write();
    h_Electron_pt->Rebin(nbins);
    h_Electron_pt->Scale(1/1000.);
    h_Electron_pt->Write();
    h_Electron_eta->Rebin(nbins);
    h_Electron_eta->Scale(1/1000.);
    h_Electron_eta->Write();
    leading_lepton_pt->Rebin(nbins);
    leading_lepton_pt->Scale(1/1000.);
    leading_lepton_pt->Write();
    ofile->Write();
}

int main(int argc, char **argv)
{
    string inputFile = argv[1];
    string outputFile = argv[2];

    TFile *ifile = TFile::Open(inputFile.c_str());
    TFile *ofile = TFile::Open(outputFile.c_str(), "RECREATE");
    TH1D *h = (TH1D*)ifile->Get("Muon_Electron_invariant_Mass");
    h->Rebin(100);
    ofile->Write();
}