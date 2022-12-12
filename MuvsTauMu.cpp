#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

// include user defined histograms and auxiliary macros
#include "Python_Analysis/corrections/roccor/RoccoR.cc"

// correctionlib
#include "correction.h"
using namespace std;
using correction::CorrectionSet;

#define MAX_ARRAY_SIZE 128
#define GEN_MAX_ARRAY_SIZE 1024

bool isFromTau(int size, Int_t *GenId, Int_t *GenParent, int initialID);

void Mixed_Analysis(string inputFile, string ofile){

cout<<"Call completed!"<<endl;

    TH1D* h_Muon_pt = new TH1D("h_Muon_pt","h_Muon_pt",100,0,200);
    TH1D* h_Muon_eta = new TH1D("h_Muon_eta","h_Muon_eta",100,-5,5);
    TH1D* h_Muon_sip = new TH1D("h_Muon_sip", "h_Muon_sip", 500, 0, 35);
    TH1D* h_Muon_ip = new TH1D("h_Muon_ip", "h_Muon_ip", 100, 0,0.02);
    TH1D* h_Muon_dxy = new TH1D("h_Muon_dxy","h_Muon_dxy",120,-0.03,0.03);
    TH1D* h_Muon_dxyE = new TH1D("h_Muon_dxyE","h_Muon_dxyE",50,0,0.01);
    TH1D* h_Muon_dxyS = new TH1D("h_Muon_dxyS","h_Muon_dxyS",60,0,20);
    TH1D* h_Muon_dz = new TH1D("h_Muon_dz","h_Muon_dz",250,-0.05,0.05);
    TH1D* h_Muon_dzE = new TH1D("h_Muon_dzE","h_Muon_dzE",100,0,0.05);    
    TH1D* h_Muon_dzS = new TH1D("h_Muon_dzS","h_Muon_dzS",200,0,20);
	// define the taumu histograms
    TH1D* h_Tau_pt = new TH1D("h_Tau_pt","h_Tau_pt",100,0,200);
    TH1D* h_Tau_eta = new TH1D("h_Tau_eta","h_Tau_eta",100,-5,5);
    TH1D* h_Tau_sip = new TH1D("h_Tau_sip", "h_Tau_sip", 500, 0, 35);
    TH1D* h_Tau_ip = new TH1D("h_Tau_ip", "h_Tau_ip", 100, 0,0.02);
    TH1D* h_Tau_dxy = new TH1D("h_Tau_dxy","h_Tau_dxy",120,-0.03,0.03);
    TH1D* h_Tau_dxyE = new TH1D("h_Tau_dxyE","h_Tau_dxyE",50,0,0.01);
    TH1D* h_Tau_dxyS = new TH1D("h_Tau_dxyS","h_Tau_dxyS",60,0,20);
    TH1D* h_Tau_dz = new TH1D("h_Tau_dz","h_Tau_dz",250,-0.05,0.05);
    TH1D* h_Tau_dzE = new TH1D("h_Tau_dzE","h_Tau_dzE",100,0,0.05);    
    TH1D* h_Tau_dzS = new TH1D("h_Tau_dzS","h_Tau_dzS",200,0,20);

   
h_Muon_pt->Sumw2(); h_Muon_eta->Sumw2(); h_Muon_sip->Sumw2(); h_Muon_ip->Sumw2(); h_Muon_dxy->Sumw2();
h_Muon_dxyE->Sumw2(); h_Muon_dxyS->Sumw2(); h_Muon_dz->Sumw2(); h_Muon_dzE->Sumw2(); h_Muon_dzS->Sumw2();
	
h_Tau_pt->Sumw2(); h_Tau_eta->Sumw2(); h_Tau_sip->Sumw2(); h_Tau_ip->Sumw2(); h_Tau_dxy->Sumw2();
h_Tau_dxyE->Sumw2(); h_Tau_dxyS->Sumw2(); h_Tau_dz->Sumw2(); h_Tau_dzE->Sumw2(); h_Tau_dzS->Sumw2();

    TFile *fin = TFile::Open(("root://cms-xrd-global.cern.ch/"+inputFile).c_str());

    TTree *tin = static_cast<TTree *>(fin->Get("Events"));

    // Set all branches to 0
    tin->SetBranchStatus("*", 0);
    // get the pt
    Float_t Muon_pt[MAX_ARRAY_SIZE], Electron_pt[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_pt", 1);
    tin->SetBranchAddress("Electron_pt", &Electron_pt);
    tin->SetBranchStatus("Muon_pt", 1);
    tin->SetBranchAddress("Muon_pt", &Muon_pt);
    // get the number of muons, electrons
    UInt_t nMuon, nElectron;
    tin->SetBranchStatus("nElectron", 1);
    tin->SetBranchAddress("nElectron", &nElectron);
    tin->SetBranchStatus("nMuon", 1);
    tin->SetBranchAddress("nMuon", &nMuon);
    // get the eta
    Float_t Muon_eta[MAX_ARRAY_SIZE], Electron_eta[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_eta", 1);
    tin->SetBranchAddress("Electron_eta", &Electron_eta);
    tin->SetBranchStatus("Muon_eta", 1);
    tin->SetBranchAddress("Muon_eta", &Muon_eta);
    // get the phi
    Float_t Muon_phi[MAX_ARRAY_SIZE], Electron_phi[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_phi", 1);
    tin->SetBranchAddress("Electron_phi", &Electron_phi);
    tin->SetBranchStatus("Muon_phi", 1);
    tin->SetBranchAddress("Muon_phi", &Muon_phi);
    // get the mass
    Float_t Muon_mass[MAX_ARRAY_SIZE], Electron_mass[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_mass", 1);
    tin->SetBranchAddress("Electron_mass", &Electron_mass);
    tin->SetBranchStatus("Muon_mass", 1);
    tin->SetBranchAddress("Muon_mass", &Muon_mass);

    // get gen quantities
    Int_t Muon_genPartIdx[MAX_ARRAY_SIZE], Electron_genPartIdx[MAX_ARRAY_SIZE];
    Int_t GenPart_pdgId[GEN_MAX_ARRAY_SIZE], GenPart_genPartIdxMother[GEN_MAX_ARRAY_SIZE], Jet_genJetIdx[MAX_ARRAY_SIZE];
    UChar_t Muon_genPartFlav[MAX_ARRAY_SIZE], Electron_genPartFlav[MAX_ARRAY_SIZE];
    UInt_t nGenPart;
    Float_t GenPart_pt[GEN_MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_genPartIdx", 1);
    tin->SetBranchStatus("Electron_genPartFlav", 1);
    tin->SetBranchStatus("Muon_genPartIdx", 1);
    tin->SetBranchStatus("Muon_genPartFlav", 1);
    tin->SetBranchStatus("GenPart_pdgId", 1);
    tin->SetBranchStatus("GenPart_genPartIdxMother", 1);
    tin->SetBranchStatus("nGenPart", 1);
    tin->SetBranchStatus("Jet_genJetIdx",1);
    tin->SetBranchStatus("GenPart_pt",1);
    tin->SetBranchAddress("nGenPart", &nGenPart);
    tin->SetBranchAddress("Electron_genPartIdx", &Electron_genPartIdx);
    tin->SetBranchAddress("Electron_genPartFlav", &Electron_genPartFlav);
    tin->SetBranchAddress("Muon_genPartIdx", &Muon_genPartIdx);
    tin->SetBranchAddress("Muon_genPartFlav", &Muon_genPartFlav);
    tin->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId);
    tin->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
    tin->SetBranchAddress("Jet_genJetIdx",&Jet_genJetIdx);
    tin->SetBranchAddress("GenPart_pt",&GenPart_pt);
    // collect the trigger information
    Bool_t HLT_IsoMu24, HLT_Ele32_WPTight_Gsf;
    tin->SetBranchStatus("HLT_IsoMu24", 1);
    tin->SetBranchStatus("HLT_Ele32_WPTight_Gsf", 1);
    tin->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);
    tin->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);

    // collect the triggger Ids
    Int_t Muon_charge[MAX_ARRAY_SIZE], Electron_charge[MAX_ARRAY_SIZE],Muon_nTrackerLayers[MAX_ARRAY_SIZE];
    Bool_t Electron_mvaFall17V2Iso_WP90[MAX_ARRAY_SIZE], Muon_triggerIdLoose[MAX_ARRAY_SIZE], Muon_tightId[MAX_ARRAY_SIZE];
    Float_t Muon_pfRelIso04_all[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_tightId", 1);
    tin->SetBranchStatus("Muon_charge", 1);
    tin->SetBranchStatus("Muon_triggerIdLoose", 1);
    tin->SetBranchStatus("Muon_pfRelIso04_all", 1);
    tin->SetBranchStatus("Electron_charge", 1);
    tin->SetBranchStatus("Electron_mvaFall17V2Iso_WP90", 1);
    tin->SetBranchStatus("Muon_nTrackerLayers", 1);
    tin->SetBranchAddress("Electron_mvaFall17V2Iso_WP90", &Electron_mvaFall17V2Iso_WP90);
    tin->SetBranchAddress("Muon_tightId", &Muon_tightId);
    tin->SetBranchAddress("Muon_charge", &Muon_charge);
    tin->SetBranchAddress("Muon_triggerIdLoose", &Muon_triggerIdLoose);
    tin->SetBranchAddress("Muon_pfRelIso04_all", &Muon_pfRelIso04_all);
    tin->SetBranchAddress("Electron_charge", &Electron_charge);
    tin->SetBranchAddress("Muon_nTrackerLayers", &Muon_nTrackerLayers);

    //mu study branches
    Float_t Muon_dxy[MAX_ARRAY_SIZE],Muon_dxyErr[MAX_ARRAY_SIZE],Muon_dz[MAX_ARRAY_SIZE],Muon_dzErr[MAX_ARRAY_SIZE];
    Float_t Muon_ip3d[MAX_ARRAY_SIZE],Muon_sip3d[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_dxy", 1);
    tin->SetBranchStatus("Muon_dxyErr", 1);
    tin->SetBranchStatus("Muon_dz", 1);
    tin->SetBranchStatus("Muon_dzErr", 1);
    tin->SetBranchStatus("Muon_ip3d", 1);
    tin->SetBranchStatus("Muon_sip3d", 1);
    tin->SetBranchAddress("Muon_dxy", &Muon_dxy);
    tin->SetBranchAddress("Muon_dxyErr", &Muon_dxyErr);
    tin->SetBranchAddress("Muon_dz", &Muon_dz);
    tin->SetBranchAddress("Muon_dzErr", &Muon_dzErr);
    tin->SetBranchAddress("Muon_ip3d", &Muon_ip3d);
    tin->SetBranchAddress("Muon_sip3d", &Muon_sip3d);

    // pu stuff
    Float_t N_pu_vertices;
    tin->SetBranchStatus("Pileup_nTrueInt", 1);
    tin->SetBranchAddress("Pileup_nTrueInt", &N_pu_vertices);

    int n_dropped = 0;
    int trigger_dropped = 0;
    UInt_t nEv = tin->GetEntries();
    unsigned int n_events = nEv;
    TLorentzVector *Muon_p4 = new TLorentzVector();
    TLorentzVector *Electron_p4 = new TLorentzVector();

    // allow pt, inv mass, and eta to be stored in a Branch
    Float_t leading_lepton_pt, invMass, electron_eta, electron_pt, muon_eta, muon_pt;
    Float_t muon_eta_from_W, muon_pt_from_W, electron_eta_from_W, electron_pt_from_W;
    float Weight;

    // open correctionfiles
    
    string muon_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/muon_Z.json.gz";
    string electron_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/electron.json.gz";
    string pileup_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/puWeights.json.gz";
    
    auto muon_c_set = CorrectionSet::from_file(muon_json);
    auto ele_c_set = CorrectionSet::from_file(electron_json);
    auto pu_c_set = CorrectionSet::from_file(pileup_json);

    auto muon_trigger = muon_c_set->at("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight");
    auto muon_id = muon_c_set->at("NUM_TightID_DEN_genTracks");
    auto muon_iso = muon_c_set->at("NUM_TightRelIso_DEN_TightIDandIPCut");
    auto electron_id = ele_c_set->at("UL-Electron-ID-SF");
    auto pu_correction = pu_c_set->at("Collisions18_UltraLegacy_goldenJSON");
    
    TFile *fecorr_trig = new TFile("/afs/cern.ch/user/g/gdamolin/public/Riccardo_egammaTriggerEfficiency_2018_20200422.root");
    TH2F * EleTrigHisto= static_cast<TH2F *>(fecorr_trig->Get("EGamma_SF2D")); 
   
    TRandom3 * RndGen=new TRandom3();

    RoccoR rc;
    rc.init("/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/roccor/RoccoR2018UL.txt");
    
    // save the histograms in a new File

    TFile *fout = new TFile(ofile.c_str(), "RECREATE");
    

    #pragma omp parallel for
    for (UInt_t i = 0; i <nEv; i++)
    {
        tin->GetEntry(i);
        if (i % 100000 == 0) std::cout << "Processing entry " << i << " of " << nEv << endl;
        
	// apply triggers
        if (!(HLT_IsoMu24 || HLT_Ele32_WPTight_Gsf)){
            trigger_dropped++;
            continue;
        };

        Int_t muon_idx = -1;
        for (UInt_t j = 0; j < nMuon; j++){
            if ((Muon_pt[j] > 27. && abs(Muon_eta[j]) < 2.4 && Muon_tightId[j] && Muon_pfRelIso04_all[j] < 0.15)){
                muon_idx = j;
		int NMCparticle=Muon_genPartIdx[muon_idx];
		double scmMC;
		if(NMCparticle>=0) {
			scmMC=rc.kSpreadMC(Muon_charge[muon_idx],Muon_pt[muon_idx],Muon_eta[muon_idx],Muon_phi[muon_idx],GenPart_pt[NMCparticle]);
			}
		else {
			scmMC=rc.kSmearMC(Muon_charge[muon_idx],Muon_pt[muon_idx],Muon_eta[muon_idx],Muon_phi[muon_idx],Muon_nTrackerLayers[muon_idx],RndGen->Rndm());
			}
		Muon_pt[muon_idx]*= scmMC;
                Muon_p4->SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], Muon_mass[j]);
		if(Muon_p4->Pt()<26) { muon_idx = -1; continue;} //if after rochester below pT threshold of trigger SF, reject muon
                else break;
            }
        }
        if (muon_idx==-1)  {
            n_dropped++;
            continue;
        }
        Weight = 1;
        Weight *= pu_correction->evaluate({N_pu_vertices, "nominal"});
			
        if(HLT_IsoMu24) {Weight *= muon_trigger->evaluate({"2018_UL", abs(Muon_eta[muon_idx]), Muon_pt[muon_idx], "sf"});} 
        Weight *= muon_id->evaluate({"2018_UL", abs(Muon_eta[muon_idx]), Muon_pt[muon_idx], "sf"}); 
        Weight *= muon_iso->evaluate({"2018_UL", abs(Muon_eta[muon_idx]), Muon_pt[muon_idx], "sf"});

	 
        Int_t electron_idx = -1;
        for (UInt_t j = 0; j < nElectron; j++){
            if ((Electron_pt[j] > 35 && abs(Electron_eta[j]) < 2.4 && Electron_mvaFall17V2Iso_WP90[j])){
		if((abs(Electron_eta[j])>1.44) && (abs(Electron_eta[j])<1.57)) {continue;} //remove electrons in the acceptance break
                Electron_p4->SetPtEtaPhiM(Electron_pt[j], Electron_eta[j], Electron_phi[j], Electron_mass[j]);
                if (Electron_p4->DeltaR(*Muon_p4) < 0.4) {continue;}
                else{electron_idx = j; break;}
            }
        }
        if (electron_idx==-1) {
            n_dropped++;
            continue;
        }

        Weight *= electron_id->evaluate({"2018", "sf", "wp90iso", abs(Electron_eta[electron_idx]), Electron_pt[electron_idx]}); 
        Weight *= electron_id->evaluate({"2018", "sf", "RecoAbove20", abs(Electron_eta[electron_idx]), Electron_pt[electron_idx]});
	if(HLT_Ele32_WPTight_Gsf) {
            //retrieve Histo
            int bin = EleTrigHisto->FindBin(Electron_eta[electron_idx],Electron_pt[electron_idx]);
            float temp= EleTrigHisto->GetBinContent(bin);
            Weight*=temp;
            }
        bool selection = ((muon_idx > -1) && (electron_idx > -1));
       
        selection = selection && (Muon_charge[muon_idx] * Electron_charge[electron_idx]) < 0;    
	
        if (!selection){ n_dropped++;  continue;}
	bool mufromtau=isFromTau(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, Muon_genPartIdx[muon_idx]);
	//bool efromtau=isFromTau(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, Electron_genPartIdx[electron_idx]);

	if (mufromtau){
		h_Tau_pt->Fill(Muon_p4->Pt(),Weight); h_Tau_eta->Fill(Muon_p4->Eta(),Weight); h_Tau_sip->Fill(Muon_sip3d[muon_idx],Weight); 
		h_Tau_ip->Fill(Muon_ip3d[muon_idx],Weight); h_Tau_dxy->Fill(Muon_dxy[muon_idx],Weight);
		h_Tau_dxyE->Fill(Muon_dxyErr[muon_idx],Weight); h_Tau_dxyS->Fill(Muon_dxy[muon_idx]/Muon_dxyErr[muon_idx],Weight);
		h_Tau_dz->Fill(Muon_dz[muon_idx],Weight); h_Tau_dzE->Fill(Muon_dzErr[muon_idx],Weight);
		h_Tau_dzS->Fill(Muon_dz[muon_idx]/Muon_dzErr[muon_idx],Weight);
	}
	else{
		h_Muon_pt->Fill(Muon_p4->Pt(),Weight); h_Muon_eta->Fill(Muon_p4->Eta(),Weight); h_Muon_sip->Fill(Muon_sip3d[muon_idx],Weight);
		h_Muon_ip->Fill(Muon_ip3d[muon_idx],Weight); h_Muon_dxy->Fill(Muon_dxy[muon_idx],Weight);
		h_Muon_dxyE->Fill(Muon_dxyErr[muon_idx],Weight); h_Muon_dxyS->Fill(Muon_dxy[muon_idx]/Muon_dxyErr[muon_idx],Weight);
		h_Muon_dz->Fill(Muon_dz[muon_idx],Weight); h_Muon_dzE->Fill(Muon_dzErr[muon_idx],Weight);
		h_Muon_dzS->Fill(Muon_dz[muon_idx]/Muon_dzErr[muon_idx],Weight);
	}
			

    }

    delete fecorr_trig;

    std::cout << "NeV = " << nEv << endl;
    std::cout << "trigger dropped = " << trigger_dropped << endl;
    std::cout << "selections dropped = " << n_dropped << endl; //remember the cross trigger in Data

    std::cout << "Fraction of events discarded by trigger = " << (trigger_dropped * 1. / nEv) << endl;
    int Rem_trigger=nEv-trigger_dropped; //remember the cross trigger in Data
    std::cout << "Fraction of events removed by selections = " << (n_dropped * 1. / Rem_trigger) << endl;
    std::cout << "Final number of events "<< Rem_trigger - n_dropped<<endl;

    h_Muon_pt->Write(); h_Muon_eta->Write(); h_Muon_sip->Write(); h_Muon_ip->Write(); h_Muon_dxy->Write();
    h_Muon_dxyE->Write(); h_Muon_dxyS->Write(); h_Muon_dz->Write(); h_Muon_dzE->Write(); h_Muon_dzS->Write();
	
    h_Tau_pt->Write(); h_Tau_eta->Write(); h_Tau_sip->Write(); h_Tau_ip->Write(); h_Tau_dxy->Write();
    h_Tau_dxyE->Write(); h_Tau_dxyS->Write(); h_Tau_dz->Write(); h_Tau_dzE->Write(); h_Tau_dzS->Write();
    
    fout->Write();
    fout->Close();
    
}

bool isFromTau(int size, Int_t *GenId, Int_t *GenParent, int initialID){
	if (initialID < 0){return false;}
	// retrieve first PDG ID number
	int startPdg = GenId[initialID];
	int newID = initialID, newPdg = startPdg;
	// look for the parent; if the parent is of same PDGID of starting particle, iterate until parent is different particle
	while (newPdg == startPdg)	{
		if (newID > size)  {std::cout << "WARNING: index " << newID << " exceeding max size " << size << std::endl;}
		newID = GenParent[newID];
		newPdg = GenId[newID];
		if (abs(newPdg) == 15)	return true;
	}
	return false;
}

int main(int argc, char **argv){

    string inputFile = argv[1];
    auto const pos = inputFile.find_last_of('/');
    string inname=inputFile.substr(pos + 1);
    string outputFile = "/afs/cern.ch/user/g/gdamolin/public/tempB/MuvsTMu/Out"+inname;

    Mixed_Analysis(inputFile, outputFile);

    return 0;
}
