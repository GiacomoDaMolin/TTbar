#ifndef Auxiliary_cpp
#define Auxiliary_cpp
#include <iostream>
// inputs are: size of arrays GenID and GenPar
// the pointer at the start of the array of GenParticles_PDGID
// the pointer at the start of the array of GenParticles_Parent
// initialID: last index of the particle (ex: Electron_genPartIdx[j]=initialID)

// should be called from func as isFromW(MAX_ARRAY_SIZE,GenPart_pdgId,GenPart_genPartIdxMother,Electron_genPartIdx[j])
//arguments: (array of PdgID of Gen Particles (Int_t), array of GenParent of Gen Particles (Int_t, contains the index to the parent of the selected GenParticle
// InitialID  is the index of the starting muon in the GenParticles array)
bool isFromW(int size, Int_t *GenId, Int_t *GenParent, int initialID)
{
	if (initialID < 0)
	{
		return false;
	}
	// retrieve first PDG ID number
	int startPdg = GenId[initialID];
	int newID = initialID, newPdg = startPdg;
	// look for the parent; if the parent is of same PDGID of starting particle, iterate until parent is different particle
	while (newPdg == startPdg)
	{
		if (newID > size)
		{
			std::cout << "WARNING: index " << newID << " exceeding max size " << size << std::endl;
		}
		newID = GenParent[newID];
		newPdg = GenId[newID];
		if (abs(newPdg) == 24)
			return true;
	}
	return false;
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

double InvertPhi(double phi){
 double invphi=phi+M_PI;
 if (invphi>M_PI){invphi=invphi-2*M_PI;}
 return invphi;
}

void HistIniz(){

 h_Muon_pt->Sumw2();
    h_Muon_eta->Sumw2();
    h_Electron_pt->Sumw2();
    h_Electron_eta->Sumw2();
    h_Muon_pt_weighted->Sumw2();
    h_Muon_eta_weighted->Sumw2();
    h_Electron_pt_weighted->Sumw2();
    h_Electron_eta_weighted->Sumw2();
    h_Muon_pt_from_W->Sumw2();
    h_Muon_eta_from_W->Sumw2();
    h_Electron_pt_from_W->Sumw2();
    h_Electron_eta_from_W->Sumw2();
    h_Muon_pt_weighted_from_W->Sumw2();
    h_Muon_eta_weighted_from_W->Sumw2();
    h_Electron_pt_weighted_from_W->Sumw2();
    h_Electron_eta_weighted_from_W->Sumw2();
    h_Muon_pt_trigger->Sumw2();
    h_Muon_eta_trigger->Sumw2();
    h_Electron_pt_trigger->Sumw2();
    h_Electron_eta_trigger->Sumw2();
    h_Muon_Electron_invariant_mass->Sumw2();
    h_Muon_Muon_invariant_mass->Sumw2();
    h_Electron_Electron_invariant_mass->Sumw2();
    h_Muon_Electron_invariant_mass_weighted->Sumw2();
    h_Muon_Muon_invariant_mass_weighted->Sumw2();
    h_Electron_Electron_invariant_mass_weighted->Sumw2();
    h_leading_lepton_pt->Sumw2();
    h_leading_lepton_pt_weighted->Sumw2();
    h_LooseJets->Sumw2();
    h_MediumJets->Sumw2();
    h_TightJets->Sumw2();
	h_dR_allJets->Sumw2();
	h_dR_lbJets->Sumw2();
	h_dR_mbJets ->Sumw2();
      h_Apl_allJets->Sumw2();
	h_Apl_lbJets->Sumw2();
	h_Apl_mbJets->Sumw2();
	h_Phi_allJets->Sumw2();
	h_Phi_lbJets->Sumw2();
	h_Phi_mbJets->Sumw2();
	h_acopla_emu->Sumw2();
	h_NJets->Sumw2();

}
#endif // Auxiliary_cpp
