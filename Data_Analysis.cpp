#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

// include user defined histograms and auxiliary macros
#include "Histodef.cpp"
#include "Auxiliary.cpp"
#include "roccor/RoccoR.cc"

using namespace std;

#define MAX_ARRAY_SIZE 128

void DataAnalysis(string inputFile, string ofile, bool IsFirstDataSet)
{

    TFile *fin = TFile::Open(inputFile.c_str());
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

    // collect the trigger information
    // triggers were 27 and 35 before
    // Triggers now changed to looser cuts due to Michele's request
    Bool_t HLT_IsoMu24, HLT_Ele32_WPTight_Gsf;
    tin->SetBranchStatus("HLT_IsoMu24", 1);
    tin->SetBranchStatus("HLT_Ele32_WPTight_Gsf", 1);
    tin->SetBranchAddress("HLT_IsoMu24", &HLT_IsoMu24);
    tin->SetBranchAddress("HLT_Ele32_WPTight_Gsf", &HLT_Ele32_WPTight_Gsf);

    // collect the triggger Ids
    Int_t Muon_charge[MAX_ARRAY_SIZE], Electron_charge[MAX_ARRAY_SIZE];
    Bool_t Electron_mvaFall17V2Iso_WP90[MAX_ARRAY_SIZE], Muon_triggerIdLoose[MAX_ARRAY_SIZE], Muon_tightId[MAX_ARRAY_SIZE];
    Float_t Muon_pfRelIso04_all[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Muon_tightId", 1);
    tin->SetBranchStatus("Muon_charge", 1);
    tin->SetBranchStatus("Muon_triggerIdLoose", 1);
    tin->SetBranchStatus("Muon_pfRelIso04_all", 1);
    tin->SetBranchStatus("Electron_charge", 1);
    tin->SetBranchStatus("Electron_mvaFall17V2Iso_WP90", 1);
    tin->SetBranchAddress("Electron_mvaFall17V2Iso_WP90", &Electron_mvaFall17V2Iso_WP90);
    tin->SetBranchAddress("Muon_tightId", &Muon_tightId);
    tin->SetBranchAddress("Muon_charge", &Muon_charge);
    tin->SetBranchAddress("Muon_triggerIdLoose", &Muon_triggerIdLoose);
    tin->SetBranchAddress("Muon_pfRelIso04_all", &Muon_pfRelIso04_all);
    tin->SetBranchAddress("Electron_charge", &Electron_charge);

    // Jet tagging , FlavB is the recomennded one, DeepB was used by Anup
    Float_t Jet_btagDeepFlavB[MAX_ARRAY_SIZE], Jet_btagDeepB[MAX_ARRAY_SIZE];
    UInt_t nJet;
    tin->SetBranchStatus("Jet_btagDeepB", 1);
    tin->SetBranchStatus("Jet_btagDeepFlavB", 1);
    tin->SetBranchStatus("nJet", 1);
    tin->SetBranchAddress("nJet", &nJet);
    tin->SetBranchAddress("Jet_btagDeepFlavB", &Jet_btagDeepFlavB);
    tin->SetBranchAddress("Jet_btagDeepB", &Jet_btagDeepB);

    int non_matching_muon = 0, non_matching_electron = 0;
    int n_dropped = 0;
    int trigger_dropped = 0,crosstrigger=0;
    const auto nEv = tin->GetEntries();
    TLorentzVector *Muon_p4 = new TLorentzVector();
    TLorentzVector *Electron_p4 = new TLorentzVector();
    TLorentzVector *MainBjet_p4 = new TLorentzVector();
    TLorentzVector *OppositeBjet_p4 = new TLorentzVector();

       // allow pt, inv mass, and eta to be stored in a Branch
    Float_t leading_lepton_pt, invMass, electron_eta, electron_pt, muon_eta, muon_pt;
    Float_t muon_eta_from_W, muon_pt_from_W, electron_eta_from_W, electron_pt_from_W;
    TFile *fout =new TFile(ofile.c_str(),"RECREATE");
    // create a new tree for the output
    TTree *tout = new TTree("tout","tout");
    // set the branches for the output tree
    tout->Branch("leading_lepton_pt", &leading_lepton_pt);
    tout->Branch("invMass", &invMass);
    tout->Branch("electron_eta", &electron_eta);
    tout->Branch("electron_pt", &electron_pt);
    tout->Branch("muon_eta", &muon_eta);
    tout->Branch("muon_pt", &muon_pt);

    int Nloose=0, Nmedium=0, Ntight=0;
    float dR_muE,dR_mujet,dR_ejet,dR_allJets,dR_lbJets,dR_mbJets,Apl_allJets,Apl_lbJets,Apl_mbJets,Phi_allJets,Phi_lbJets,Phi_mbJets, PTbjet;

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

    RoccoR rc;
    rc.init("/afs/cern.ch/user/g/gdamolin/Johan/TTbar/RoccoR2018UL.txt");
    
    for (UInt_t i = 0; i < nEv; i++)
    {
        tin->GetEntry(i);
        if (i % 100000 == 0)
            std::cout << "Processing entry " << i << " of " << nEv << std::endl;
        // apply triggers

        if (!(HLT_IsoMu24 || HLT_Ele32_WPTight_Gsf))
        {
            trigger_dropped++;
            continue;
        };

        // avoid cross triggers
        if (!IsFirstDataSet && HLT_Ele32_WPTight_Gsf && HLT_IsoMu24)
        {
            crosstrigger++;
            continue;
        }

        // loop over the muons and electrons and only keep the fist ones that pass the requirements
        Int_t muon_idx = -1;
        for (UInt_t j = 0; j < nMuon; j++){
            if ((Muon_pt[j]>27. && abs(Muon_eta[j])<2.4 && Muon_tightId[j] && Muon_pfRelIso04_all[j] < 0.15)){
                muon_idx = j;
                Muon_p4->SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], Muon_mass[j]);
                break;
            }
        }
        Int_t electron_idx = -1;
        for (UInt_t j = 0; j < nElectron; j++){
            if ((Electron_pt[j]>35 && abs(Electron_eta[j])<2.4 && Electron_mvaFall17V2Iso_WP90[j])){
                Electron_p4->SetPtEtaPhiM(Electron_pt[j], Electron_eta[j], Electron_phi[j], Electron_mass[j]);
		if(Electron_p4->DeltaR(*Muon_p4)<0.4) {continue;}
                else {electron_idx = j; break;}
            }
        }
        bool selection = ((muon_idx > -1) && (electron_idx > -1));
        // check the seleected objects for opposite charge
        selection = selection && (Muon_charge[muon_idx] * Electron_charge[electron_idx]) < 0;
        // the tight working point is 0.71, medium 0.2783, loose 0.0490
        Float_t jet_btag_deepFlav_wp = 0.2783;
        // the wp are: (0.1355, 0.4506, 0.7738)
        // Float_t jet_btag_deep_wp = 0.4506;
        // cycle through btags and check if one passes the tagging WP
        bool one_Bjet = false;
        int id_m_jet=-1;
	Nloose=0, Nmedium=0, Ntight=0;
        for (size_t j = 0; j < nJet; j++){
            if (Jet_btagDeepFlavB[j] > jet_btag_deepFlav_wp){
		if(!one_Bjet) {
			MainBjet_p4->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
			OppositeBjet_p4->SetPtEtaPhiM(Jet_pt[j], -1*Jet_eta[j], InvertPhi(Jet_phi[j]), Jet_mass[j]);
				
			if ( MainBjet_p4->DeltaR(*Muon_p4)<0.4 || MainBjet_p4->DeltaR(*Electron_p4)<0.4) continue;
                	one_Bjet = true; id_m_jet=j;
                	}	
            }
            if (Jet_btagDeepFlavB[j] > 0.0490)	Nloose++;
	    if (Jet_btagDeepFlavB[j] > 0.2783)	Nmedium++;
	    if (Jet_btagDeepFlavB[j] > 0.71)	Ntight++;
        }
        selection = selection && (one_Bjet);
        if (!selection)
        {
            n_dropped++;
            continue;
        }
        PTbjet=MainBjet_p4->Pt();

	double scmDT=rc.kScaleDT(Muon_charge[muon_idx],Muon_pt[muon_idx],Muon_eta[muon_idx],Muon_phi[muon_idx]);
        Muon_pt[muon_idx]*= scmDT;
	Muon_p4->SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], Muon_mass[j]);
        

        dR_mujet=Muon_p4->DeltaR(*MainBjet_p4);
	dR_ejet=Electron_p4->DeltaR(*MainBjet_p4);
	dR_muE=Muon_p4->DeltaR(*Electron_p4);

        // check whether muon or electron is the leading one
        if (Muon_p4->Pt() > Electron_p4->Pt()){
            // fill the hist
            leading_lepton_pt = Muon_p4->Pt();  
        } else {
            leading_lepton_pt = Electron_p4->Pt();
            
        }
        h_leading_lepton_pt->Fill(leading_lepton_pt);
        // fill the histograms
        muon_pt = Muon_pt[muon_idx];
        muon_eta = Muon_eta[muon_idx];
        electron_pt = Electron_pt[electron_idx];
        electron_eta = Electron_eta[electron_idx];

        h_Muon_pt->Fill(muon_pt);
        h_Muon_eta->Fill(muon_eta);

        h_Electron_pt->Fill(electron_pt);
        h_Electron_eta->Fill(electron_eta);

        
	dR_allJets=999, dR_lbJets=999, dR_mbJets=999;
	Apl_allJets=1.1,Apl_lbJets=1.1,Apl_mbJets=1.1;
	bool ok1=false,ok2=false ,ok3=false;
	for (size_t j = 0; j < nJet; j++){
		if (j==id_m_jet) continue;
		TLorentzVector *tempJet = new TLorentzVector();
		tempJet->SetPtEtaPhiM(Jet_pt[j], Jet_eta[j], Jet_phi[j], Jet_mass[j]);
		double temp=OppositeBjet_p4->DeltaR(*tempJet);

		TVector3 A(tempJet->X(),tempJet->Y(),tempJet->Z());	
		TVector3 B(MainBjet_p4->X(),MainBjet_p4->Y(),MainBjet_p4->Z());

		double tempApl=A.Dot(B)/(A.Mag()*B.Mag());

		if(temp<dR_allJets) {dR_allJets=temp; ok1=true;}
		if(tempApl<Apl_allJets) {Apl_allJets=tempApl;}
		if(Jet_btagDeepFlavB[j] > 0.0490){ ok2=true;
			if (temp<dR_lbJets){dR_lbJets=temp;}
			if (tempApl<Apl_lbJets) {Apl_lbJets=tempApl;}
			}
		if(Jet_btagDeepFlavB[j] > 0.2783){ ok3=true;
			if (temp<dR_mbJets){dR_mbJets=temp;}
			if (tempApl<Apl_mbJets) {Apl_mbJets=tempApl;}
			}
  
			
		delete tempJet;	
        	}

//dphi
	Phi_allJets=999, Phi_lbJets=999, Phi_mbJets=999;
	for (size_t j = 0; j < nJet; j++){
			if(j==id_m_jet) continue;
			double temp=Jet_phi[j]-OppositeBjet_p4->Phi();
			if (temp<-1*M_PI) temp+=2*M_PI;
			if (temp>M_PI) temp-=2*M_PI;
			if (temp<0) temp*=(-1);

			if(temp<Phi_allJets) {Phi_allJets=temp;}
			if((Jet_btagDeepFlavB[j] > 0.0490) && (temp<Phi_lbJets)) {Phi_lbJets=temp;}
			if((Jet_btagDeepFlavB[j] > 0.2783) && (temp<Phi_mbJets)) {Phi_mbJets=temp;}
        		}

        if (muon_idx > -1 && electron_idx > -1)
        {
           // calculate the invariant mass of the two muons
            invMass = (*(Muon_p4) + *(Electron_p4)).M();
            // fill the invariant mass histogram
            h_Muon_Electron_invariant_mass->Fill(invMass);
        }
	tout->Fill();
    }
    std::cout << "Total number of events: " << nEv << std::endl;
    std::cout << "Number of events discarded by trigger = " << trigger_dropped << std::endl;
    std:cout <<"Number of event removed by cross-triggering filter"<<crosstrigger<<std::endl;
    trigger_dropped+=crosstrigger;
    std::cout << "Number of events discarded by selection = " << n_dropped << std::endl;

    std::cout << "Number of events passing triggers = " << (nEv - trigger_dropped) << std::endl;
    std::cout << "Number of events passing selection = " << (nEv - trigger_dropped) - n_dropped << std::endl;

    std::cout << "Selected events over triggered events = " << (nEv - n_dropped) * 1. / (nEv - trigger_dropped) << std::endl;
    // Write the histograms to the file
    h_Muon_eta->Write();
    h_Muon_pt->Write();

    h_Electron_eta->Write();
    h_Electron_pt->Write();

    h_Muon_Electron_invariant_mass->Write();
    h_leading_lepton_pt->Write();

    fout->Write();
    fout->Close();
}

int main(int argc, char **argv)
{
    string inputFile = argv[1];
    string outputFile = argv[2];
    string boolstr=argv[3];
    bool IsFirstDataset= (boolstr=="true");
    if (IsFirstDataset) {std::cout<<"############## It is first dataset! ##################"<<std::endl;}
    if (!IsFirstDataset) {std::cout<<"@@@@@@@@@@@@@@ NOT first dataset! @@@@@@@@@@@@@@@@@@"<<std::endl;}
    DataAnalysis(inputFile, outputFile, IsFirstDataset);
}
