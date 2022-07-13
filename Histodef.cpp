#include "TH1.h"

 // Define the histogram objects
    TH1D* h_Muon_pt = new TH1D("Muon_pt", "Muon_pt", 100, 0, 200);
    TH1D* h_Muon_eta = new TH1D("Muon_eta", "Muon_eta", 100, -5, 5);
    TH1D* h_Electron_pt = new TH1D("Electron_pt", "Electron_pt", 100, 0, 200);
    TH1D* h_Electron_eta = new TH1D("Electron_eta", "Electron_eta", 100, -5, 5);

    TH1D* h_Muon_pt_trigger = new TH1D("Muon_pt_trigger", "Muon_pt_trigger", 100, 0, 200);
    TH1D* h_Muon_eta_trigger = new TH1D("Muon_eta_trigger", "Muon_eta_trigger", 100, -5, 5);
    TH1D* h_Electron_pt_trigger = new TH1D("Electron_pt_trigger", "Electron_pt_trigger", 100, 0, 200);
    TH1D* h_Electron_eta_trigger = new TH1D("Electron_eta_trigger", "Electron_eta_trigger", 100, -5, 5);


    TH1D* h_Muon_Electron_invariant_mass = new TH1D("Muon_Electron_invariant_mass", "Muon_Electron_invariant_mass", 100, 0, 200);
    TH1D* h_Muon_Muon_invariant_mass = new TH1D("Muon_Muon_invariant_mass", "Muon_Muon_invariant_mass", 100, 0, 200);
    TH1D* h_Electron_Electron_invariant_mass = new TH1D("Electron_Electron_invariant_mass", "Electron_Electron_invariant_mass", 100, 0, 200);
