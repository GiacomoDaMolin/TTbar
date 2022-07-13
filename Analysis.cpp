#include <iostream>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TLorentzVector.h"

//include user defined histograms and auxiliary macros
#include "Auxiliary.cpp"
#include "Histodef.cpp"

using namespace std;

#define MAX_ARRAY_SIZE 128


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
    Int_t Muon_genPartIdx[MAX_ARRAY_SIZE], Electron_genPartIdx[MAX_ARRAY_SIZE], GenPart_pdgId[MAX_ARRAY_SIZE*4], GenPart_genPartIdxMother[MAX_ARRAY_SIZE*4],; //These last guys are actually huge, better be careful!
    UChar_t Muon_genPartFlav[MAX_ARRAY_SIZE], Electron_genPartFlav[MAX_ARRAY_SIZE];
    tin->SetBranchStatus("Electron_genPartIdx", 1);
    tin->SetBranchStatus("Electron_genPartFlav", 1);
    tin->SetBranchStatus("Muon_genPartIdx", 1);
    tin->SetBranchStatus("Muon_genPartFlav", 1);
    tin->SetBranchStatus("GenPart_pdgId", 1);
    tin->SetBranchStatus("GenPart_genPartIdxMother",1);
    tin->SetBranchAddress("Electron_genPartIdx", &Electron_genPartIdx);
    tin->SetBranchAddress("Electron_genPartFlav", &Electron_genPartFlav);
    tin->SetBranchAddress("Muon_genPartIdx", &Muon_genPartIdx);
    tin->SetBranchAddress("Muon_genPartFlav", &Muon_genPartFlav);
    tin->SetBranchAddress("GenPart_pdgId", &GenPart_pdgId);
    tin->SetBranchAddress("GenPart_genPartIdxMother",&GenPart_genPartIdxMother);

    // collect the trigger information
    Bool_t HLT_IsoMu27, HLT_Ele35_WPTight_Gsf; 
    tin->SetBranchStatus("HLT_IsoMu27", 1);
    tin->SetBranchStatus("HLT_Ele35_WPTight_Gsf", 1);
    tin->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27);
    tin->SetBranchAddress("HLT_Ele35_WPTight_Gsf", &HLT_Ele35_WPTight_Gsf);

    //tin->Draw("Muon_pt >> h_Muon_pt");
    /*tin->Draw("Muon_eta >> h_Muon_eta","", "GOFF");
    tin->Draw("Electron_pt >> h_Electron_pt","", "GOFF");
    tin->Draw("Electron_eta >> h_Electron_eta","", "GOFF");*/
    //tin->Draw("Muon_pt >> h_Muon_pt_trigger","(HLT_IsoMu27) || (HLT_Ele35_WPTight_Gsf)", "GOFF");
    //tin->Draw("Muon_eta >> h_Muon_eta_trigger", "(HLT_IsoMu27) || (HLT_Ele35_WPTight_Gsf)","GOFF");
    //tin->Draw("Electron_pt >> h_Electron_pt_trigger", "(HLT_IsoMu27) || (HLT_Ele35_WPTight_Gsf)","GOFF");
    //tin->Draw("Electron_eta >> h_Electron_eta_trigger", "(HLT_IsoMu27) || (HLT_Ele35_WPTight_Gsf)","GOFF");
    const auto nEv = tin->GetEntries();
    for (size_t i = 0; i < 100; i++){
        tin->GetEntry(i);
        // apply triggers
        //if (HLT_IsoMu27==0&&HLT_Ele35_WPTight_Gsf==0){
        //    continue;
        //}
        // initialise exactly 2 LorentzVectors for e and muon as there should not be more
        vector <TLorentzVector> Muon_p4, Electron_p4;
        for (int j = 0; j<nMuon; j++){
            h_Muon_pt->Fill(Muon_pt[j]);
            h_Muon_eta->Fill(Muon_eta[j]);
            TLorentzVector Muon_p4_temp;
            if (HLT_IsoMu27 || HLT_Ele35_WPTight_Gsf){
                h_Muon_pt_trigger->Fill(Muon_pt[j]);
                h_Muon_eta_trigger->Fill(Muon_eta[j]);
            }
            // match the muon to the PID of the W boson (PID=24)
            if (GenPart_pdgId[Muon_genPartIdx[j]] == 24 || GenPart_pdgId[Muon_genPartIdx[j]] == -24){ //isFromW(MAX_ARRAY_SIZE,GenPart_pdgId,GenPart_genPartIdxMother,Muon_genPartIdx[j])
                
                Muon_p4_temp.SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], Muon_mass[j]);
                Muon_p4.push_back(Muon_p4_temp);
            }
            // match the electron to the PID of the W boson (PID=24)
        }
        
        for (int j = 0; j<nElectron; j++){
            h_Electron_pt->Fill(Electron_pt[j]);
            h_Electron_eta->Fill(Electron_eta[j]);
            TLorentzVector Electron_p4_temp;

            if (HLT_IsoMu27 || HLT_Ele35_WPTight_Gsf){
                h_Electron_pt_trigger->Fill(Electron_pt[j]);
                h_Electron_eta_trigger->Fill(Electron_eta[j]);
            }
            if (GenPart_pdgId[Electron_genPartIdx[j]] == 24 || GenPart_pdgId[Electron_genPartIdx[j]] == -24){ //isFromW(MAX_ARRAY_SIZE,GenPart_pdgId,GenPart_genPartIdxMother,Electron_genPartIdx[j])
                
                Electron_p4_temp.SetPtEtaPhiM(Electron_pt[j], Electron_eta[j], Electron_phi[j], Electron_mass[j]);
                Electron_p4.push_back(Electron_p4_temp);
            }
        bool selection=false;
        //selection= (Muon_pt[0]>SOMETHING && abs(Muon_eta[0])<2.4) && (idem with ele)
        //selection = selection && opposite charge mu and e
        //selection = selection && ((muon is trigger) || (ele is trigger))
        //selection = selection && (one b-jet)
        //if (selection) {fill histo}
        }
        // check the number of muons and electrons
        cout << "Muon_p4.size() = " << Muon_p4.size() << endl;
        cout << "Electron_p4.size() = " << Electron_p4.size() << endl;
        if (Muon_p4.size() == 1 || Electron_p4.size() == 1){
            // calculate the invariant mass of the two muons
            float_t lepton_invariant_mass = (Muon_p4[0] + Electron_p4[0]).M();
            // fill the invariant mass histogram
            h_Muon_Electron_invariant_mass->Fill(lepton_invariant_mass);
        }

        if (Muon_p4.size() == 2){
            // calculate the invariant mass of the two muons
            float_t lepton_invariant_mass = (Muon_p4[0] + Muon_p4[1]).M();
            // fill the invariant mass histogram
            h_Muon_Muon_invariant_mass->Fill(lepton_invariant_mass);
        }
        if (Electron_p4.size() == 2){
            // calculate the invariant mass of the two muons
            float_t lepton_invariant_mass = (Electron_p4[0] + Electron_p4[1]).M();
            // fill the invariant mass histogram
            h_Electron_Electron_invariant_mass->Fill(lepton_invariant_mass);
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

    h_Muon_Electron_invariant_mass->Write();
    h_Muon_Muon_invariant_mass->Write();
    h_Electron_Electron_invariant_mass->Write();

    fout->Write();
    fout->Close();
    
}
