#include <string>

#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "THStack.h"
void Normalize(TH1D *a, int nbins = 10)
{
	a->Scale(1. / a->Integral());
	a->Rebin(nbins);
}

void Plotter()
{

	// retrieve histograms
	TFile *FileSig = new TFile("Out_Signal.root");
	TH1D *MuPT_sig = static_cast<TH1D *>(FileSig->Get("h_Muon_pt"));
	TH1D *MuETA_sig = static_cast<TH1D *>(FileSig->Get("h_Muon_eta"));
	TH1D *ElePT_sig = static_cast<TH1D *>(FileSig->Get("h_Electron_pt"));
	TH1D *EleETA_sig = static_cast<TH1D *>(FileSig->Get("h_Electron_eta"));
	TH1D *DiLeptMass_sig = static_cast<TH1D *>(FileSig->Get("Muon_Electron_invariant_mass"));
	Normalize(MuPT_sig);
	Normalize(MuETA_sig);
	Normalize(ElePT_sig);
	Normalize(EleETA_sig);
	Normalize(DiLeptMass_sig);

	TFile *FileDY = new TFile("Out_DY_all.root");
	TH1D *MuPT_DY = static_cast<TH1D *>(FileDY->Get("h_Muon_pt"));
	TH1D *MuETA_DY = static_cast<TH1D *>(FileDY->Get("h_Muon_eta"));
	TH1D *ElePT_DY = static_cast<TH1D *>(FileDY->Get("h_Electron_pt"));
	TH1D *EleETA_DY = static_cast<TH1D *>(FileDY->Get("h_Electron_eta"));
	TH1D *DiLeptMass_DY = static_cast<TH1D *>(FileDY->Get("Muon_Electron_invariant_mass"));
	Normalize(MuPT_DY);
	Normalize(MuETA_DY);
	Normalize(ElePT_DY);
	Normalize(EleETA_DY);
	Normalize(DiLeptMass_DY);

	TFile *FileW = new TFile("Out_W_all.root");
	TH1D *MuPT_W = static_cast<TH1D *>(FileW->Get("h_Muon_pt"));
	TH1D *MuETA_W = static_cast<TH1D *>(FileW->Get("h_Muon_eta"));
	TH1D *ElePT_W = static_cast<TH1D *>(FileW->Get("h_Electron_pt"));
	TH1D *EleETA_W = static_cast<TH1D *>(FileW->Get("h_Electron_eta"));
	TH1D *DiLeptMass_W = static_cast<TH1D *>(FileW->Get("Muon_Electron_invariant_mass"));
	Normalize(MuPT_W);
	Normalize(MuETA_W);
	Normalize(ElePT_W);
	Normalize(EleETA_W);
	Normalize(DiLeptMass_W);

	TFile *FileSt = new TFile("Out_St_all.root");
	TH1D *MuPT_St = static_cast<TH1D *>(FileSt->Get("h_Muon_pt"));
	TH1D *MuETA_St = static_cast<TH1D *>(FileSt->Get("h_Muon_eta"));
	TH1D *ElePT_St = static_cast<TH1D *>(FileSt->Get("h_Electron_pt"));
	TH1D *EleETA_St = static_cast<TH1D *>(FileSt->Get("h_Electron_eta"));
	TH1D *DiLeptMass_St = static_cast<TH1D *>(FileSt->Get("Muon_Electron_invariant_mass"));
	Normalize(MuPT_St);
	Normalize(MuETA_St);
	Normalize(ElePT_St);
	Normalize(EleETA_St);
	Normalize(DiLeptMass_St);

	std::cout << "Red = Signal;   Blue=DY;   Green=W+jets;    Black=Single top" << std::endl;
	// plot histograms
	TCanvas *d1 = new TCanvas("MuPT");
	THStack *hs1 = new THStack("hs1", "");
	MuPT_sig->SetXTitle("p_{T} of the Muon [GeV]");
	MuPT_sig->SetLineColor(kRed);
	MuPT_DY->SetLineColor(kBlue);
	MuPT_W->SetLineColor(kGreen);
	MuPT_St->SetLineColor(kBlack);
	MuPT_sig->SetLineWidth(3);
	MuPT_DY->SetLineWidth(3);
	MuPT_W->SetLineWidth(3);
	MuPT_St->SetLineWidth(3);
	hs1->Add(MuPT_sig);
	hs1->Add(MuPT_DY);
	hs1->Add(MuPT_W);
	hs1->Add(MuPT_St);
	hs1->Draw("nostack");
	gPad->BuildLegend(0.7, 0.8, 1., 1.);
	gPad->Update();

	TCanvas *d2 = new TCanvas("MuETA");
	THStack *hs2 = new THStack("hs2", "");
	MuETA_sig->SetXTitle("Eta of the Muon");
	MuETA_sig->SetLineColor(kRed);
	MuETA_DY->SetLineColor(kBlue);
	MuETA_W->SetLineColor(kGreen);
	MuETA_St->SetLineColor(kBlack);
	MuETA_sig->SetLineWidth(3);
	MuETA_DY->SetLineWidth(3);
	MuETA_W->SetLineWidth(3);
	MuETA_St->SetLineWidth(3);
	hs2->Add(MuETA_sig);
	hs2->Add(MuETA_DY);
	hs2->Add(MuETA_W);
	hs2->Add(MuETA_St);
	hs2->Draw("nostack");
	gPad->BuildLegend(0.7, 0.8, 1., 1.);
	gPad->Update();

	TCanvas *d3 = new TCanvas("ElePT");
	THStack *hs3 = new THStack("hs3", "");
	ElePT_sig->SetXTitle("p_{T} of the Electron [GeV]");
	ElePT_sig->SetLineColor(kRed);
	ElePT_DY->SetLineColor(kBlue);
	ElePT_W->SetLineColor(kGreen);
	ElePT_St->SetLineColor(kBlack);
	ElePT_sig->SetLineWidth(3);
	ElePT_DY->SetLineWidth(3);
	ElePT_W->SetLineWidth(3);
	ElePT_St->SetLineWidth(3);
	hs3->Add(ElePT_sig);
	hs3->Add(ElePT_DY);
	hs3->Add(ElePT_W);
	hs3->Add(ElePT_St);
	hs3->Draw("nostack");
	gPad->BuildLegend(0.7, 0.8, 1., 1.);
	gPad->Update();

	TCanvas *d4 = new TCanvas("EleETA");
	THStack *hs4 = new THStack("hs4", "");
	EleETA_sig->SetXTitle("eta of the Electron");
	EleETA_sig->SetLineColor(kRed);
	EleETA_DY->SetLineColor(kBlue);
	EleETA_W->SetLineColor(kGreen);
	EleETA_St->SetLineColor(kBlack);
	EleETA_sig->SetLineWidth(3);
	EleETA_DY->SetLineWidth(3);
	EleETA_W->SetLineWidth(3);
	EleETA_St->SetLineWidth(3);
	hs4->Add(EleETA_sig);
	hs4->Add(EleETA_DY);
	hs4->Add(EleETA_W);
	hs4->Add(EleETA_St);
	hs4->Draw("nostack");
	gPad->BuildLegend(0.7, 0.8, 1., 1.);
	gPad->Update();

	TCanvas *d5 = new TCanvas("Mass");
	THStack *hs5 = new THStack("hs5", "");
	DiLeptMass_sig->SetXTitle("Mass of the e-mu system [GeV]");
	DiLeptMass_sig->SetLineColor(kRed);
	DiLeptMass_DY->SetLineColor(kBlue);
	DiLeptMass_W->SetLineColor(kGreen);
	DiLeptMass_St->SetLineColor(kBlack);
	DiLeptMass_sig->SetLineWidth(3);
	DiLeptMass_DY->SetLineWidth(3);
	DiLeptMass_W->SetLineWidth(3);
	DiLeptMass_St->SetLineWidth(3);
	hs5->Add(DiLeptMass_sig);
	hs5->Add(DiLeptMass_DY);
	hs5->Add(DiLeptMass_W);
	hs5->Add(DiLeptMass_St);
	hs5->Draw("nostack");
	gPad->BuildLegend(0.7, 0.8, 1., 1.);
	gPad->Update();

	return;
}
