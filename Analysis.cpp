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
#define GEN_MAX_ARRAY_SIZE 1024


void Analysis(string inputFile, string ofile);
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
    Analysis(inputFile, outputFile);
}

void Analysis(string inputFile, string ofile){

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
    Int_t Muon_genPartIdx[MAX_ARRAY_SIZE], Electron_genPartIdx[MAX_ARRAY_SIZE];
    Int_t GenPart_pdgId[GEN_MAX_ARRAY_SIZE], GenPart_genPartIdxMother[GEN_MAX_ARRAY_SIZE]; //These last guys are actually huge, better be careful!
    UChar_t Muon_genPartFlav[MAX_ARRAY_SIZE], Electron_genPartFlav[MAX_ARRAY_SIZE];
    UInt_t nGenPart;
    tin->SetBranchStatus("Electron_genPartIdx", 1);
    tin->SetBranchStatus("Electron_genPartFlav", 1);
    tin->SetBranchStatus("Muon_genPartIdx", 1);
    tin->SetBranchStatus("Muon_genPartFlav", 1);
    tin->SetBranchStatus("GenPart_pdgId", 1);
    tin->SetBranchStatus("GenPart_genPartIdxMother",1);
    tin->SetBranchStatus("nGenPart", 1);
    tin->SetBranchAddress("nGenPart", &nGenPart);
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
    vector <TLorentzVector*> Muon_p4(20, nullptr), Electron_p4(20, nullptr);
    for (size_t i = 0; i < nEv; i++){
        //std::cout << "Event = " << i << std::endl;
        tin->GetEntry(i);
        // apply triggers
        //if (HLT_IsoMu27==0&&HLT_Ele35_WPTight_Gsf==0){
        //    continue;
        //}
        // initialise exactly 2 LorentzVectors for e and muon as there should not be more
        size_t nMuon_p4 = 0, nElectron_p4 = 0;
        for (int j = 0; j<nMuon; j++){
            h_Muon_pt->Fill(Muon_pt[j]);
            h_Muon_eta->Fill(Muon_eta[j]);
            if (HLT_IsoMu27 || HLT_Ele35_WPTight_Gsf){
                h_Muon_pt_trigger->Fill(Muon_pt[j]);
                h_Muon_eta_trigger->Fill(Muon_eta[j]);
            }
            // match the muon to the PID of the W boson (PID=24)
            //printMCTree(nGenPart, GenPart_pdgId,GenPart_genPartIdxMother, Muon_genPartIdx[j]);
            if (isFromW(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother,Muon_genPartIdx[j])){
                if (Muon_p4[nMuon_p4] == nullptr){
                    Muon_p4[nMuon_p4] = new TLorentzVector();
                }
                // nMuon_p4++ increments by one and returns the previuos value
                //std::cout << "Muon " << nMuon_p4 << std::endl;
                Muon_p4[nMuon_p4++]->SetPtEtaPhiM(Muon_pt[j], Muon_eta[j], Muon_phi[j], Muon_mass[j]);
            }
        }
        
        for (int j = 0; j<nElectron; j++){
            h_Electron_pt->Fill(Electron_pt[j]);
            h_Electron_eta->Fill(Electron_eta[j]);

            if (HLT_IsoMu27 || HLT_Ele35_WPTight_Gsf){
                h_Electron_pt_trigger->Fill(Electron_pt[j]);
                h_Electron_eta_trigger->Fill(Electron_eta[j]);
            }
            //printMCTree(nGenPart, GenPart_pdgId,GenPart_genPartIdxMother, Electron_genPartIdx[j]);
            if (isFromW(nGenPart,GenPart_pdgId,GenPart_genPartIdxMother,Electron_genPartIdx[j])){                 
                if (Electron_p4[nElectron_p4] == nullptr){
                    Electron_p4[nElectron_p4] = new TLorentzVector();
                }
                //std::cout << "Electron " << nElectron_p4 << std::endl;
                Electron_p4[nElectron_p4++]->SetPtEtaPhiM(Electron_pt[j], Electron_eta[j], Electron_phi[j], Electron_mass[j]);
            }
            bool selection=false;
            //selection= (Muon_pt[0]>SOMETHING && abs(Muon_eta[0])<2.4) && (idem with ele)
            //selection = selection && opposite charge mu and e
            //selection = selection && ((muon is trigger) || (ele is trigger))
            //selection = selection && (one b-jet)
            //if (selection) {fill histo}
        }
        // check the number of muons and electrons
        //cout << "Muon_p4.size() = " << nMuon_p4 << endl;
        //cout << "Electron_p4.size() = " << nElectron_p4 << endl;
        if (nMuon_p4 == 1 && nElectron_p4 == 1){
            // calculate the invariant mass of the two muons
            float_t lepton_invariant_mass = (*(Muon_p4[0]) + *(Electron_p4[0])).M();
            // fill the invariant mass histogram
            h_Muon_Electron_invariant_mass->Fill(lepton_invariant_mass);
        }

        if (nMuon_p4 == 2){
            // calculate the invariant mass of the two muons
            float_t lepton_invariant_mass = (*(Muon_p4[0]) + *(Muon_p4[1])).M();
            // fill the invariant mass histogram
            h_Muon_Muon_invariant_mass->Fill(lepton_invariant_mass);
        }
        if (nElectron_p4 == 2){
            // calculate the invariant mass of the two muons
            float_t lepton_invariant_mass = (*(Electron_p4[0]) + *(Electron_p4[1])).M();
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
