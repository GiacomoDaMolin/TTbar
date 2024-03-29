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
#include "Histodef.cpp"
#include "Auxiliary.cpp"
#include "Python_Analysis/corrections/roccor/RoccoR.cc"

// correctionlib
#include "correction.h"
using namespace std;
using correction::CorrectionSet;

#define MAX_ARRAY_SIZE 128
#define GEN_MAX_ARRAY_SIZE 1024

// function to calculate the weight for each event
// the weight is calculated as the product of luminosity and cross section of the process times the genWeight,
// LATER TO BE divided by the number of generated events OF ALL FILES OF THE DATASET(S)
double getWeight(double luminosity, double crossSection, Float_t genWeight, double SumWeights)
{
    return (luminosity * crossSection * genWeight); // / SumWeights;
}

void Mixed_Analysis(string inputFile, string ofile, double crossSection = -1, double IntLuminosity = 59.827879506, bool Signal = false)
{

    if (crossSection < 0. || IntLuminosity < 0.)
    {
        std::cout << "WARNING: crossection " << crossSection << " and Integrated luminosity " << IntLuminosity << endl;
    }

cout<<"Call completed!"<<endl;

    TFile *fin = TFile::Open(inputFile.c_str());
    TTree *trun = static_cast<TTree *>(fin->Get("Runs"));
    Long64_t genEventCount;
    Double_t genEventSumw;
    trun->SetBranchStatus("*", 0);
    trun->SetBranchStatus("genEventSumw", 1);
    trun->SetBranchStatus("genEventCount", 1);
    trun->SetBranchAddress("genEventSumw", &genEventSumw);
    trun->SetBranchAddress("genEventCount", &genEventCount);

    trun->GetEntry(0);

    TTree *tin = static_cast<TTree *>(fin->Get("Events"));

    // Set all branches to 0
    tin->SetBranchStatus("*", 0);
    // get the pt
    Float_t Muon_pt[MAX_ARRAY_SIZE], Electron_pt[MAX_ARRAY_SIZE], Jet_pt[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_pt", 1);
    tin->SetBranchAddress("Electron_pt", &Electron_pt);
    tin->SetBranchStatus("Muon_pt", 1);
    tin->SetBranchAddress("Muon_pt", &Muon_pt);
    tin->SetBranchStatus("Jet_pt", 1);
    tin->SetBranchAddress("Jet_pt", &Jet_pt);
    // get the number of muons, electrons
    UInt_t nMuon, nElectron;
    tin->SetBranchStatus("nElectron", 1);
    tin->SetBranchAddress("nElectron", &nElectron);
    tin->SetBranchStatus("nMuon", 1);
    tin->SetBranchAddress("nMuon", &nMuon);
    // get the eta
    Float_t Muon_eta[MAX_ARRAY_SIZE], Electron_eta[MAX_ARRAY_SIZE], Jet_eta[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_eta", 1);
    tin->SetBranchAddress("Electron_eta", &Electron_eta);
    tin->SetBranchStatus("Muon_eta", 1);
    tin->SetBranchAddress("Muon_eta", &Muon_eta);
    tin->SetBranchStatus("Jet_eta", 1);
    tin->SetBranchAddress("Jet_eta", &Jet_eta);
    // get the phi
    Float_t Muon_phi[MAX_ARRAY_SIZE], Electron_phi[MAX_ARRAY_SIZE], Jet_phi[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_phi", 1);
    tin->SetBranchAddress("Electron_phi", &Electron_phi);
    tin->SetBranchStatus("Muon_phi", 1);
    tin->SetBranchAddress("Muon_phi", &Muon_phi);
    tin->SetBranchStatus("Jet_phi", 1);
    tin->SetBranchAddress("Jet_phi", &Jet_phi);
    // get the mass
    Float_t Muon_mass[MAX_ARRAY_SIZE], Electron_mass[MAX_ARRAY_SIZE], Jet_mass[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_mass", 1);
    tin->SetBranchAddress("Electron_mass", &Electron_mass);
    tin->SetBranchStatus("Muon_mass", 1);
    tin->SetBranchAddress("Muon_mass", &Muon_mass);
    tin->SetBranchStatus("Jet_mass", 1);
    tin->SetBranchAddress("Jet_mass", &Jet_mass);

    // get gen quantities
    Int_t Muon_genPartIdx[MAX_ARRAY_SIZE], Electron_genPartIdx[MAX_ARRAY_SIZE];
    Int_t GenPart_pdgId[GEN_MAX_ARRAY_SIZE], GenPart_genPartIdxMother[GEN_MAX_ARRAY_SIZE], GenPart_statusFlags[GEN_MAX_ARRAY_SIZE], Jet_genJetIdx[MAX_ARRAY_SIZE];
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
    tin->SetBranchStatus("GenPart_statusFlags",1);
    tin->SetBranchAddress("nGenPart", &nGenPart);
    tin->SetBranchAddress("Electron_genPartIdx", &Electron_genPartIdx);
    tin->SetBranchAddress("Electron_genPartFlav", &Electron_genPartFlav);
    tin->SetBranchAddress("Muon_genPartIdx", &Muon_genPartIdx);
    tin->SetBranchAddress("Muon_genPartFlav", &Muon_genPartFlav);
    tin->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId);
    tin->SetBranchAddress("GenPart_genPartIdxMother", &GenPart_genPartIdxMother);
    tin->SetBranchAddress("Jet_genJetIdx",&Jet_genJetIdx);
    tin->SetBranchAddress("GenPart_pt",&GenPart_pt);
    tin->SetBranchAddress("GenPart_statusFlags",&GenPart_statusFlags);

    // collect the trigger information
    Bool_t HLT_IsoMu24, HLT_Ele32_WPTight_Gsf;
    tin->SetBranchStatus("HLT_IsoMu24", 1);
    tin->SetBranchStatus("HLT_Ele32_WPTight_Gsf", 1);
    tin->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);
    tin->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);

    // collect the triggger Ids
    Int_t Muon_charge[MAX_ARRAY_SIZE], Electron_charge[MAX_ARRAY_SIZE],Muon_nTrackerLayers[MAX_ARRAY_SIZE];
    Bool_t Electron_mvaFall17V2Iso_WP90[MAX_ARRAY_SIZE], Muon_triggerIdLoose[MAX_ARRAY_SIZE], Muon_tightId[MAX_ARRAY_SIZE];
    Float_t Muon_pfRelIso04_all[MAX_ARRAY_SIZE], Electron_ip3d[MAX_ARRAY_SIZE], Electron_sip3d[MAX_ARRAY_SIZE], Electron_dxy[MAX_ARRAY_SIZE], Electron_dz[MAX_ARRAY_SIZE], Muon_dxy[MAX_ARRAY_SIZE], Muon_dz[MAX_ARRAY_SIZE], Muon_ip3d[MAX_ARRAY_SIZE], Muon_sip3d[MAX_ARRAY_SIZE];

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

    tin->SetBranchStatus("Electron_ip3d", 1);
    tin->SetBranchStatus("Electron_sip3d", 1);
    tin->SetBranchStatus("Electron_dxy", 1);
    tin->SetBranchStatus("Electron_dz", 1);
    tin->SetBranchStatus("Muon_dxy", 1);
    tin->SetBranchStatus("Muon_dz", 1);
    tin->SetBranchStatus("Muon_ip3d", 1);
    tin->SetBranchStatus("Muon_sip3d", 1);
    tin->SetBranchAddress("Electron_dz", &Electron_dz);
    tin->SetBranchAddress("Electron_ip3d", &Electron_ip3d);
    tin->SetBranchAddress("Electron_sip3d", &Electron_sip3d);
    tin->SetBranchAddress("Electron_dxy", &Electron_dxy);
    tin->SetBranchAddress("Muon_dxy", &Muon_dxy);
    tin->SetBranchAddress("Muon_dz", &Muon_dz);
    tin->SetBranchAddress("Muon_ip3d", &Muon_ip3d);
    tin->SetBranchAddress("Muon_sip3d", &Muon_sip3d);

    // Jet tagging and ID, FlavB is the recomended one, DeepB was used by Anup
    Float_t Jet_btagDeepFlavB[MAX_ARRAY_SIZE], Jet_btagDeepB[MAX_ARRAY_SIZE];
    UInt_t nJet;
    Int_t Jet_jetId[MAX_ARRAY_SIZE], Jet_puId[MAX_ARRAY_SIZE],Jet_hadronFlavour[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Jet_btagDeepB", 1);
    tin->SetBranchStatus("Jet_btagDeepFlavB", 1);
    tin->SetBranchStatus("nJet", 1);
    tin->SetBranchStatus("Jet_jetId", 1);
    tin->SetBranchStatus("Jet_puId", 1);
    tin->SetBranchStatus("Jet_hadronFlavour", 1);
    tin->SetBranchAddress("nJet", &nJet);
    tin->SetBranchAddress("Jet_btagDeepFlavB", &Jet_btagDeepFlavB);
    tin->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB);   
    tin->SetBranchAddress("Jet_jetId", &Jet_jetId);
    tin->SetBranchAddress("Jet_puId", &Jet_puId);
    tin->SetBranchAddress("Jet_hadronFlavour", &Jet_hadronFlavour);

    // pu stuff
    Float_t N_pu_vertices;
    tin->SetBranchStatus("Pileup_nTrueInt", 1);
    tin->SetBranchAddress("Pileup_nTrueInt", &N_pu_vertices);

    // gen weight
    Float_t genWeight;
    tin->SetBranchStatus("genWeight", 1);
    tin->SetBranchAddress("genWeight", &genWeight);


    //L1
    Float_t L1PreFiringWeight_Nom;
    tin->SetBranchStatus("L1PreFiringWeight_Nom", 1);
    tin->SetBranchAddress("L1PreFiringWeight_Nom", &L1PreFiringWeight_Nom);
    Int_t TrigObj_id[MAX_ARRAY_SIZE]; UInt_t nTrigObj;
    Float_t TrigObj_eta[MAX_ARRAY_SIZE], TrigObj_phi[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("TrigObj_id", 1);
    tin->SetBranchStatus("TrigObj_eta", 1);
    tin->SetBranchStatus("TrigObj_phi", 1);
    tin->SetBranchStatus("nTrigObj", 1);
    tin->SetBranchAddress("TrigObj_id", &TrigObj_id);
    tin->SetBranchAddress("TrigObj_eta", &TrigObj_eta);
    tin->SetBranchAddress("TrigObj_phi", &TrigObj_phi); 
    tin->SetBranchAddress("nTrigObj", &nTrigObj);
    int evenottrigMatch=0;

    int non_matching_muon = 0, non_matching_electron = 0;
    int n_dropped = 0;
    int trigger_dropped = 0;
    UInt_t nEv = tin->GetEntries();
    unsigned int n_events = nEv;
    TLorentzVector *Muon_p4 = new TLorentzVector();
    TLorentzVector *Electron_p4 = new TLorentzVector();
    TLorentzVector *MainBjet_p4 = new TLorentzVector();
    TLorentzVector *OppositeBjet_p4 = new TLorentzVector();

    // allow pt, inv mass, and eta to be stored in a Branch
    Float_t leading_lepton_pt, invMass, electron_eta, electron_pt, muon_eta, muon_pt;
    Float_t muon_eta_from_W, muon_pt_from_W, electron_eta_from_W, electron_pt_from_W;
    float Weight;

    // open correctionfiles
    
    string muon_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/muon_Z.json.gz";
    string electron_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/electron.json.gz";
    string jets_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/jet_jmar.json";
    string b_tag_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/btagging.json.gz";
    string pileup_json = "/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/puWeights.json.gz";

    
    auto muon_c_set = CorrectionSet::from_file(muon_json);
    auto ele_c_set = CorrectionSet::from_file(electron_json);
    auto jet_c_set = CorrectionSet::from_file(jets_json);
    auto btag_c_set = CorrectionSet::from_file(b_tag_json);
    auto pu_c_set = CorrectionSet::from_file(pileup_json);

    auto muon_trigger = muon_c_set->at("NUM_IsoMu24_DEN_CutBasedIdTight_and_PFIsoTight");
    auto muon_id = muon_c_set->at("NUM_TightID_DEN_genTracks");
    auto muon_iso = muon_c_set->at("NUM_TightRelIso_DEN_TightIDandIPCut");
    auto electron_id = ele_c_set->at("UL-Electron-ID-SF");
    auto jet_pu = jet_c_set->at("PUJetID_eff");
    auto b_tag = btag_c_set->at("deepJet_mujets");
    auto b_mistag= btag_c_set->at("deepJet_incl"); //only for light jets
    auto pu_correction = pu_c_set->at("Collisions18_UltraLegacy_goldenJSON");
    
    TFile *fecorr_trig = new TFile("/afs/cern.ch/user/g/gdamolin/public/Riccardo_egammaTriggerEfficiency_2018_20200422.root");
    TH2F * EleTrigHisto= static_cast<TH2F *>(fecorr_trig->Get("EGamma_SF2D"));

    TFile *fb_eff = new TFile("/afs/cern.ch/user/g/gdamolin/public/Beff_puLoose.root");
    TH2D * l_eff= static_cast<TH2D *>(fb_eff->Get("l_jets_tagged")); 
    TH2D * c_eff= static_cast<TH2D *>(fb_eff->Get("c_jets_tagged")); 
    TH2D * b_eff= static_cast<TH2D *>(fb_eff->Get("b_jets_tagged")); 
   
    TRandom3 * RndGen=new TRandom3();
    RoccoR rc;
    rc.init("/afs/cern.ch/user/g/gdamolin/Johan/TTbar/Python_Analysis/corrections/roccor/RoccoR2018UL.txt");
     
    // save the histograms in a new File

    TFile *fout = new TFile(ofile.c_str(), "RECREATE");
    
    // create a new tree for the output
    TTree *tout = new TTree("tout", "tout");
    TTree *trun_out = new TTree("Run_out", "Run_out");
    // set the branches for the output tree
    tout->Branch("leading_lepton_pt", &leading_lepton_pt);
    tout->Branch("invMass", &invMass);
    tout->Branch("electron_eta", &electron_eta);
    tout->Branch("electron_pt", &electron_pt);
    tout->Branch("muon_eta", &muon_eta);
    tout->Branch("muon_pt", &muon_pt);
    tout->Branch("Weight", &Weight);

    int Nloose = 0, Nmedium = 0, Ntight = 0, JetsNotB=0;
    float dR_muE, dR_mujet, dR_ejet, dR_allJets, dR_lbJets, dR_mbJets, Apl_allJets, Apl_lbJets, Apl_mbJets, Phi_allJets, Phi_lbJets, Phi_mbJets, PTbjet,Acopl_emu;

    bool From2Taus=false, FromTau=false;

    tout->Branch("dR_mue", &dR_muE);
    tout->Branch("dR_mujet", &dR_mujet);
    tout->Branch("dR_ejet", &dR_ejet);
    tout->Branch("dR_allJets", &dR_allJets);
    tout->Branch("dR_lbJets", &dR_lbJets);
    tout->Branch("dR_mbJets", &dR_mbJets);
    tout->Branch("Apl_lbJets", &Apl_lbJets);
    tout->Branch("Apl_allJets", &Apl_allJets);
    tout->Branch("Apl_mbJets", &Apl_mbJets);
    tout->Branch("Phi_allJets", &Phi_allJets);
    tout->Branch("Phi_lbJets", &Phi_lbJets);
    tout->Branch("Phi_mbJets", &Phi_mbJets);
    tout->Branch("PTbjet", &PTbjet);
    tout->Branch("Nloose", &Nloose);
    tout->Branch("Nmedium", &Nmedium);
    tout->Branch("Ntight", &Ntight);
    tout->Branch("JetNotB", &JetsNotB);
    tout->Branch("Acopl_emu", &Acopl_emu);
    tout->Branch("FromTau", &FromTau);
    tout->Branch("From2Taus", &From2Taus);
    tout->Branch("MuTrigger", &HLT_IsoMu24);
    tout->Branch("ETrigger", &HLT_Ele32_WPTight_Gsf);

    trun_out->Branch("genEventSumw", &genEventSumw);
    trun_out->Branch("IntLumi", &IntLuminosity);
    trun_out->Branch("xs", &crossSection);
    trun_out->Branch("nEv", &n_events);

    trun_out->Fill(); // we already called trun->GetEntry(0);

    size_t found = ofile.find_last_of("/");
    string oname=ofile.substr(found+1);
    string path=ofile.substr(0,found);
    string Tauname=path+"/Tau_"+oname;
    TFile *foutT = new TFile(Tauname.c_str(), "RECREATE");
    TTree *toutT = new TTree("toutT", "toutT");
    TTree *trun_outT = new TTree("Run_outT", "Run_outT");
    if (Signal) {
	toutT->Branch("leading_lepton_pt", &leading_lepton_pt);
    	toutT->Branch("invMass", &invMass);
    	toutT->Branch("electron_eta", &electron_eta);
   	toutT->Branch("electron_pt", &electron_pt);
   	toutT->Branch("muon_eta", &muon_eta);
    	toutT->Branch("muon_pt", &muon_pt);
    	toutT->Branch("Weight", &Weight);
	toutT->Branch("FromTau", &FromTau);
	toutT->Branch("From2Taus", &From2Taus);
	toutT->Branch("dR_mue", &dR_muE);
	toutT->Branch("dR_mujet", &dR_mujet);
	toutT->Branch("dR_ejet", &dR_ejet);
	    toutT->Branch("dR_allJets", &dR_allJets);
	    toutT->Branch("dR_lbJets", &dR_lbJets);
	    toutT->Branch("dR_mbJets", &dR_mbJets);
	    toutT->Branch("Apl_lbJets", &Apl_lbJets);
	    toutT->Branch("Apl_allJets", &Apl_allJets);
	    toutT->Branch("Apl_mbJets", &Apl_mbJets);
	    toutT->Branch("Phi_allJets", &Phi_allJets);
	    toutT->Branch("Phi_lbJets", &Phi_lbJets);
	    toutT->Branch("Phi_mbJets", &Phi_mbJets);
	    toutT->Branch("PTbjet", &PTbjet);
	    toutT->Branch("Nloose", &Nloose);
	    toutT->Branch("Nmedium", &Nmedium);
	    toutT->Branch("Ntight", &Ntight);
	    toutT->Branch("JetNotB", &JetsNotB);
	    toutT->Branch("Acopl_emu", &Acopl_emu);
            toutT->Branch("MuTrigger", &HLT_IsoMu24);
    toutT->Branch("ETrigger", &HLT_Ele32_WPTight_Gsf);


	trun_outT->Branch("genEventSumw", &genEventSumw);
        trun_outT->Branch("IntLumi", &IntLuminosity);
        trun_outT->Branch("xs", &crossSection);
        trun_outT->Branch("nEv", &n_events);
	trun_outT->Fill();
	}
    fout->cd();

    #pragma omp parallel for
    for (UInt_t i = 0; i <nEv; i++)
    {
        tin->GetEntry(i);
        if (i % 100000 == 0) std::cout << "Processing entry " << i << " of " << nEv << endl;
        
	// apply triggers
        if (!(HLT_Ele32_WPTight_Gsf)){
            trigger_dropped++;
            continue;
        };

        Int_t muon_idx = -1;
        for (UInt_t j = 0; j < nMuon; j++)
        { if ((abs(Muon_eta[j]) < 2.4 && Muon_tightId[j] && Muon_pfRelIso04_all[j] < 0.15)){
		  int NMCparticle=Muon_genPartIdx[j];
		  double scmMC;
		  if(NMCparticle>=0) {scmMC=rc.kSpreadMC(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j],GenPart_pt[NMCparticle]);}
		  else {scmMC=rc.kSmearMC(Muon_charge[j],Muon_pt[j],Muon_eta[j],Muon_phi[j],Muon_nTrackerLayers[j],RndGen->Rndm());}

		  Muon_pt[j]*= scmMC;
		  if ( Muon_pt[j] > 27.){
		        muon_idx = j;
		        Muon_p4->SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], Muon_mass[j]);
			break;
		    }
          }
        }
        if (muon_idx==-1)  {
            n_dropped++;
            continue;
        }
        Weight = getWeight(IntLuminosity, crossSection, genWeight, genEventSumw);
	Weight *=  getTopPtWeight(GenPart_pdgId,GenPart_statusFlags,GenPart_pt,nGenPart);
	double Weight2=Weight;
        Weight *= pu_correction->evaluate({N_pu_vertices, "nominal"});
	Weight*= L1PreFiringWeight_Nom;
			
        //if(HLT_IsoMu24) {Weight *= muon_trigger->evaluate({"2018_UL", abs(Muon_eta[muon_idx]), Muon_pt[muon_idx], "sf"});} 
        Weight *= muon_id->evaluate({"2018_UL", abs(Muon_eta[muon_idx]), Muon_pt[muon_idx], "sf"}); 
        Weight *= muon_iso->evaluate({"2018_UL", abs(Muon_eta[muon_idx]), Muon_pt[muon_idx], "sf"});

	 
        Int_t electron_idx = -1;
        for (UInt_t j = 0; j < nElectron; j++){
            if ((Electron_pt[j] > 35 && abs(Electron_eta[j]) < 2.4 && Electron_mvaFall17V2Iso_WP90[j] && abs(Electron_dxy[j])<0.2 && abs(Electron_dz[j])<0.5)){
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
/*
	bool saved=false;
	for(int j=0; j<nTrigObj;j++){
		if(TrigObj_id[j]==11){
		 TLorentzVector *T_p4 = new TLorentzVector();
	         T_p4->SetPtEtaPhiM(10, TrigObj_eta[j], TrigObj_phi[j], 10);
		 if(T_p4->DeltaR(*Electron_p4)<0.1)	{saved=true; delete T_p4;
					 break;}
		 delete T_p4;
		}
		if(TrigObj_id[j]==13){
			 TLorentzVector *T_p4 = new TLorentzVector();
			 T_p4->SetPtEtaPhiM(10, TrigObj_eta[j], TrigObj_phi[j], 10);
			 if(T_p4->DeltaR(*Muon_p4)<0.1)	{saved=true; delete T_p4;
					 break;}
			 delete T_p4; 
			}
	}

	if (saved==false) {evenottrigMatch++; n_dropped++; continue;}	*/
	
	if(Signal){
		bool firstistau=isFromTau(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, Muon_genPartIdx[muon_idx]);
		bool secondistau=isFromTau(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, Electron_genPartIdx[electron_idx]);
		if (firstistau||secondistau){FromTau=true;}
		else {FromTau=false;}
		if (firstistau&&secondistau){From2Taus=true;}
		else {From2Taus=false;}
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
        // check the seleected objects for opposite charge
        selection = selection && (Muon_charge[muon_idx] * Electron_charge[electron_idx]) < 0;
        // the tight working point is 0.71, medium 0.2783, loose 0.0490
        Float_t jet_btag_deepFlav_wp = 0.2783;
        bool one_Bjet = false;
        int id_m_jet = -1;
	int njets=0;
        Nloose = 0, Nmedium = 0, Ntight = 0, JetsNotB=0;
	//vectors for applying b-tag corrections
	vector<int> njet_in_collection;
	vector<int> flavor;
	vector<bool> tagged;
	double t_weight=1.;
        for (size_t j = 0; j < nJet; j++)
        {
            if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6)){
	    TLorentzVector *Tjet_p4 = new TLorentzVector();
	    Tjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	    if((Tjet_p4->DeltaR(*Muon_p4)<0.4) || (Tjet_p4->DeltaR(*Electron_p4)<0.4)) {delete Tjet_p4; continue;}
	    else {delete Tjet_p4;}
	
            //correction for pileupID
            int MC_pu = Jet_genJetIdx[j];
            float tempSF=1.,tempEff;
            //if is pileUpjet
            if (MC_pu<0 ) {
		tempSF=1.;
            	tempEff= 0;
            	}
            //if is truly a jet
            else { if (Jet_pt[j]<=50){
            	     tempSF= jet_pu->evaluate({Jet_eta[j],Jet_pt[j],"nom", "L"});
            	     tempEff= jet_pu->evaluate({Jet_eta[j],Jet_pt[j],"MCEff", "L"});
		    }
            	}
            bool passesPUID=(Jet_puId[j]>=4);
		
            if(!(Jet_pt[j]>50 || passesPUID ))	{t_weight*=(1-tempSF*tempEff)/(1-tempEff); }
            if((Jet_pt[j]>50 || passesPUID)) { 
             if(Jet_pt[j]<=50) t_weight*=tempSF; //else you are in pT>50 case: apply no sf
              //correction for b-tag
              njet_in_collection.push_back(j);
              flavor.push_back(abs(Jet_hadronFlavour[j]));
              tagged.push_back((Jet_btagDeepFlavB[j] > jet_btag_deepFlav_wp));
	      njets++;
        
	      if (Jet_btagDeepFlavB[j] < 0.0490) JetsNotB++;
	      if (Jet_btagDeepFlavB[j] > 0.0490) Nloose++;
              if (Jet_btagDeepFlavB[j] > 0.2783) Nmedium++;
              if (Jet_btagDeepFlavB[j] > 0.71) Ntight++;
            
              if (Jet_btagDeepFlavB[j] > jet_btag_deepFlav_wp){
                 if (!one_Bjet){
                      MainBjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
                      OppositeBjet_p4->SetPtEtaPhiM(Jet_pt[j], -1 * Jet_eta[j], InvertPhi(Jet_phi[j]), Jet_mass[j]);

                      if (MainBjet_p4->DeltaR(*Muon_p4) > 0.4 && MainBjet_p4->DeltaR(*Electron_p4) > 0.4){
                      	one_Bjet = true;
                      	id_m_jet = j;
			 }
                      }
               }
            }//end if(jetpt>50 !!puid==7)
          }//end kinematic if
        }//end for
          //corrections of jets already applied 
	
        Weight*=t_weight; 
            
	for(int jj=0;jj<flavor.size();jj++){
		int convflav=flavor[jj];
		if (flavor[jj]<4) convflav==0;
		if (!(convflav==0 || convflav==4 || convflav==5)) {cout<<"Something weird in the flavor of jet"<<endl;}
		if(tagged[jj]){
			if (convflav!=0) {Weight *= b_tag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});}
			else  {Weight *= b_mistag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});}
			continue;}

		//if not tagged
		if(!tagged[jj]) {
			double Eff=1.;
			double SF=1;
			if (convflav!=0) SF=b_tag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
			else SF=b_mistag->evaluate({"central", "M", convflav, abs(Jet_eta[njet_in_collection[jj]]), Jet_pt[njet_in_collection[jj]]});
			
			//Get Eff
			if(convflav==0) {
				int bin =l_eff->FindBin(Jet_pt[njet_in_collection[jj]],abs(Jet_eta[njet_in_collection[jj]]));
				Eff=l_eff->GetBinContent(bin);
				}
			if(convflav==4) {
				int bin =c_eff->FindBin(Jet_pt[njet_in_collection[jj]],abs(Jet_eta[njet_in_collection[jj]]));
				Eff=c_eff->GetBinContent(bin);
				}
			if(convflav==5) {
				int bin =b_eff->FindBin(Jet_pt[njet_in_collection[jj]],abs(Jet_eta[njet_in_collection[jj]]));
				Eff=b_eff->GetBinContent(bin);
				}
			Weight*=(1-SF*Eff)/(1-Eff);
			}
		
		}

        selection = selection && (one_Bjet);
        if (!selection){ n_dropped++;  continue;}
        h_LooseJets->Fill(Nloose, Weight);
        h_MediumJets->Fill(Nmedium, Weight);
        h_TightJets->Fill(Ntight, Weight);
        Acopl_emu=M_PI-(Electron_p4->DeltaPhi(*Muon_p4));
        h_acopla_emu->Fill(Acopl_emu,Weight);

        PTbjet = MainBjet_p4->Pt();

        dR_mujet = Muon_p4->DeltaR(*MainBjet_p4);
        dR_ejet = Electron_p4->DeltaR(*MainBjet_p4);
        dR_muE = Muon_p4->DeltaR(*Electron_p4);


        if (Muon_p4->Pt() > Electron_p4->Pt()){
            leading_lepton_pt = Muon_p4->Pt();
            h_leading_lepton_pt->Fill(leading_lepton_pt,Weight2);
            h_leading_lepton_pt_weighted->Fill(leading_lepton_pt, Weight);
        }
        else{
            leading_lepton_pt = Electron_p4->Pt();
            h_leading_lepton_pt->Fill(leading_lepton_pt,Weight2);
            h_leading_lepton_pt_weighted->Fill(leading_lepton_pt, Weight);
        }

        // fill the histograms
        muon_pt = Muon_p4->Pt();
        muon_eta = Muon_eta[muon_idx];
        electron_pt = Electron_pt[electron_idx];
        electron_eta = Electron_eta[electron_idx];

        h_Muon_pt->Fill(muon_pt, Weight2);
        h_Muon_eta->Fill(muon_eta, Weight2);

        h_Electron_pt->Fill(electron_pt,Weight2);
        h_Electron_eta->Fill(electron_eta, Weight2);
        // fill the weighted histograms
        h_Muon_pt_weighted->Fill(muon_pt, Weight);
        h_Muon_eta_weighted->Fill(muon_eta, Weight);
        h_Electron_pt_weighted->Fill(electron_pt, Weight);
        h_Electron_eta_weighted->Fill(electron_eta, Weight);

	h_NJets->Fill(njets,Weight);
        // only for signal
        if (Signal)
        {
            // cross check which index the objects have that actually originate from the W
            size_t nMuon_p4 = 0, nElectron_p4 = 0;
            for (UInt_t j = 0; j < nMuon; j++)
            {
                // match the muon to the PID of the W boson (PID=24)
                // printMCTree(nGenPart, GenPart_pdgId,GenPart_genPartIdxMother, Muon_genPartIdx[j]);
                if (isFromW(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, Muon_genPartIdx[j]))
                {
                    muon_pt_from_W = Muon_p4->Pt();
                    muon_eta_from_W = Muon_eta[j];
                    h_Muon_pt_from_W->Fill(muon_pt_from_W,Weight2);
                    h_Muon_eta_from_W->Fill(muon_eta_from_W,Weight2);
                    h_Muon_pt_weighted_from_W->Fill(muon_pt_from_W, Weight);
                    h_Muon_eta_weighted_from_W->Fill(muon_eta_from_W, Weight);
                    if (muon_idx != j)
                        non_matching_muon++;
                }
            }

            for (UInt_t j = 0; j < nElectron; j++)
            {
                if (isFromW(nGenPart, GenPart_pdgId, GenPart_genPartIdxMother, Electron_genPartIdx[j]))
                {
                    electron_pt_from_W = Electron_pt[j];
                    electron_eta_from_W = Electron_eta[j];
                    h_Electron_pt_from_W->Fill(electron_pt_from_W,Weight2);
                    h_Electron_eta_from_W->Fill(electron_eta_from_W,Weight2);
                    h_Electron_pt_weighted_from_W->Fill(electron_pt_from_W, Weight);
                    h_Electron_eta_weighted_from_W->Fill(electron_eta_from_W, Weight);
                    if (electron_idx != j)
                        non_matching_electron++;
                }
            }
        }
        // END only for signal

        dR_allJets = 999, dR_lbJets = 999, dR_mbJets = 999;
        Apl_allJets = 1.1, Apl_lbJets = 1.1, Apl_mbJets = 1.1;
        for (size_t j = 0; j < nJet; j++)
        {  
          if (j == id_m_jet)
                continue;
	  TLorentzVector *Tjet_p4 = new TLorentzVector();
	  Tjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	  if((Tjet_p4->DeltaR(*Muon_p4)<0.4) || (Tjet_p4->DeltaR(*Electron_p4)<0.4)) {delete Tjet_p4; continue;}
	  
          if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6) && (Jet_pt[j]>50 || (Jet_puId[j]>=4))){
            TLorentzVector *tempJet = Tjet_p4;
            double temp = OppositeBjet_p4->DeltaR(*tempJet);

            TVector3 A(tempJet->X(), tempJet->Y(), tempJet->Z());
            TVector3 B(MainBjet_p4->X(), MainBjet_p4->Y(), MainBjet_p4->Z());

            double tempApl = A.Dot(B) / (A.Mag() * B.Mag());

            if (temp < dR_allJets){dR_allJets = temp;}
            if (tempApl < Apl_allJets){Apl_allJets = tempApl;}
            if (Jet_btagDeepFlavB[j] > 0.0490){
                if (temp < dR_lbJets){dR_lbJets = temp;}
                if (tempApl < Apl_lbJets){Apl_lbJets = tempApl;}
            }
            if (Jet_btagDeepFlavB[j] > 0.2783){
                if (temp < dR_mbJets){dR_mbJets = temp;}
                if (tempApl < Apl_mbJets){Apl_mbJets = tempApl;}
            }
            delete Tjet_p4;
         }//end if
        } //end for

        // dphi
        Phi_allJets = 999, Phi_lbJets = 999, Phi_mbJets = 999;
        for (size_t j = 0; j < nJet; j++){ 
	    if (j == id_m_jet) continue;
	    TLorentzVector *Tjet_p4 = new TLorentzVector();
	    Tjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
	    if((Tjet_p4->DeltaR(*Muon_p4)<0.4) || (Tjet_p4->DeltaR(*Electron_p4)<0.4)) {delete Tjet_p4; continue;}
	    else {delete Tjet_p4;}

	    if((abs(Jet_eta[j]) < 2.4) && Jet_pt[j]>25 && (Jet_jetId[j]==2 || Jet_jetId[j]==6) && (Jet_pt[j]>50 || (Jet_puId[j]>=4))){
		    double temp = Jet_phi[j] - OppositeBjet_p4->Phi();
		    if (temp < -1 * M_PI) temp += 2 * M_PI;
		    if (temp > M_PI) temp -= 2 * M_PI;
		    if (temp < 0) temp *= (-1);

		    if (temp < Phi_allJets){Phi_allJets = temp;}
		    if ((Jet_btagDeepFlavB[j] > 0.0490) && (temp < Phi_lbJets)){Phi_lbJets = temp;}
		    if ((Jet_btagDeepFlavB[j] > 0.2783) && (temp < Phi_mbJets)){Phi_mbJets = temp;}
	    }
        }

        if (muon_idx > -1 && electron_idx > -1){
            invMass = (*(Muon_p4) + *(Electron_p4)).M();
            // fill the invariant mass histogram
            h_Muon_Electron_invariant_mass->Fill(invMass, Weight2);
            h_Muon_Electron_invariant_mass_weighted->Fill(invMass, Weight);
        }
	h_Trigger->Fill(HLT_IsoMu24+2*HLT_Ele32_WPTight_Gsf,Weight); 
	h_mu_3dsig->Fill(Muon_sip3d[muon_idx],Weight);
	h_mu_3d->Fill(Muon_ip3d[muon_idx],Weight);
	h_mu_dxy->Fill(abs(Muon_dxy[muon_idx]),Weight);
	h_e_3dsig->Fill(Electron_sip3d[electron_idx],Weight);
	h_e_3d->Fill(Electron_ip3d[electron_idx],Weight);
	h_e_dxy->Fill(abs(Electron_dxy[electron_idx]),Weight);

	b_pt->Fill(MainBjet_p4->Pt(),Weight);
	jethole->Fill(MainBjet_p4->Eta(),MainBjet_p4->Phi(),Weight);
	ehole->Fill(Electron_p4->Eta(),Electron_p4->Phi(),Weight);


        // fill the tree
        if(Signal && FromTau) {toutT->Fill();}
	else {tout->Fill();}
    }

    delete fecorr_trig;
    cout<<" Is this signal?"<<endl;
    if (Signal) cout<<"Yes"<<endl;
    else cout<<"No"<<endl;
    std::cout << "non_matching_muon = " << non_matching_muon << endl;
    std::cout << "non_matching_electron = " << non_matching_electron << endl;

    cout<<"L1 requirement removes "<< evenottrigMatch << " events"<<endl;

    std::cout << "NeV = " << nEv << endl;
    std::cout << "trigger dropped = " << trigger_dropped << endl;
    std::cout << "selections dropped = " << n_dropped << endl; //remember the cross trigger in Data

    std::cout << "Fraction of events discarded by trigger = " << (trigger_dropped * 1. / nEv) << endl;
    int Rem_trigger=nEv-trigger_dropped; //remember the cross trigger in Data
    std::cout << "Fraction of events removed by selections = " << (n_dropped * 1. / Rem_trigger) << endl;
    std::cout << "Final number of events "<< Rem_trigger - n_dropped<<endl;

    tout->Write();
    trun_out->Write();
	
    // Write the histograms to the file

    HistWrite();

    fout->Close();
    if (Signal) {
	cout<<"Saving Tau File!"<<endl;
	foutT->cd();
	toutT->Write();
	trun_outT->Write();
	foutT->Write();
	foutT->Close();
	}
   else cout<<"Not writing tree in Tau File"<<endl;
}

int main(int argc, char **argv)
{

    string inputFile = argv[1];
    string outputFile = argv[2];
    double crossSection = atof(argv[3]);
    double IntLuminosity = atof(argv[4]);
    string boolstr = argv[5];
    bool Signal = (boolstr == "true" || boolstr == "True");

    HistIniz();

    Mixed_Analysis(inputFile, outputFile, crossSection, IntLuminosity, Signal);

    return 0;
}
