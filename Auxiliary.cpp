#ifndef Auxiliary_cpp
#define Auxiliary_cpp
#include <iostream>
#include<vector>
using std::vector;
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


/**
     @short computes the weight to apply to the ttbar simulation to correct for the top pT spectrum
     unless two "last copy top quarks" are found the event the weight returned is the trivial one
     cf. https://twiki.cern.ch/twiki/bin/view/CMS/TopPtReweighting
   */
  float getTopPtWeight(Int_t * pdgId,Int_t *statusFlags,Float_t * pt, Int_t Ngen) {

    float wgt=1.0;

    //lastcopy is signalled by the 13th bit: (statusFlags >> 13 ) & 0x1 == (statusFlags/8192)%2
    //auto topquarks_lastcopy = ROOT::VecOps::Nonzero( abs(pdgId)==6 && ((statusFlags/8192)%2) );
    vector<Int_t> topquarks_lastcopy;
    for(Int_t i=0;i<Ngen;i++){
	if(abs(pdgId[i])==6 && ((statusFlags[i]/8192)%2))
		topquarks_lastcopy.push_back(i);
	}
    if(topquarks_lastcopy.size()!=2) return wgt;

    //NNLO / NLO SF parameterization
    for(auto i : topquarks_lastcopy) {
      wgt *= 0.103 * exp(-0.0118 * pt[i]) - 0.000134 * pt[i] + 0.973;
    }
    return sqrt(wgt);
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
  h_vsPTandEta_onlye->Sumw2();
  h_vsPTandEta_onlymu->Sumw2();
  h_vsPTandEta_muande->Sumw2();

h_Trigger->Sumw2();
h_mu_3dsig->Sumw2();
h_mu_3d ->Sumw2();
h_mu_dxy->Sumw2();
h_e_3dsig->Sumw2();
h_e_3d->Sumw2();
h_e_dxy->Sumw2();

b_pt->Sumw2();
jethole->Sumw2();
ehole->Sumw2();

}

void HistWrite(){
		h_Muon_pt_weighted->SetBinContent(h_Muon_pt_weighted->GetNbinsX(), h_Muon_pt_weighted->GetBinContent(h_Muon_pt_weighted->GetNbinsX()) + h_Muon_pt_weighted->GetBinContent(h_Muon_pt_weighted->GetNbinsX() + 1));
	h_Electron_pt_weighted->SetBinContent(h_Electron_pt_weighted->GetNbinsX(), h_Electron_pt_weighted->GetBinContent(h_Electron_pt_weighted->GetNbinsX()) + h_Electron_pt_weighted->GetBinContent(h_Electron_pt_weighted->GetNbinsX() + 1));
	 h_Muon_Muon_invariant_mass_weighted->SetBinContent( h_Muon_Muon_invariant_mass_weighted->GetNbinsX(),  h_Muon_Muon_invariant_mass_weighted->GetBinContent( h_Muon_Muon_invariant_mass_weighted->GetNbinsX()) +  h_Muon_Muon_invariant_mass_weighted->GetBinContent( h_Muon_Muon_invariant_mass_weighted->GetNbinsX() + 1));

h_Muon_pt->Write();
    h_Muon_eta->Write();
    h_Electron_pt->Write();
    h_Electron_eta->Write();
    h_Muon_pt_weighted->Write();
    h_Muon_eta_weighted->Write();
    h_Electron_pt_weighted->Write();
    h_Electron_eta_weighted->Write();
    h_Muon_pt_from_W->Write();
    h_Muon_eta_from_W->Write();
    h_Electron_pt_from_W->Write();
    h_Electron_eta_from_W->Write();
    h_Muon_pt_weighted_from_W->Write();
    h_Muon_eta_weighted_from_W->Write();
    h_Electron_pt_weighted_from_W->Write();
    h_Electron_eta_weighted_from_W->Write();
    h_Muon_pt_trigger->Write();
    h_Muon_eta_trigger->Write();
    h_Electron_pt_trigger->Write();
    h_Electron_eta_trigger->Write();
    h_Muon_Electron_invariant_mass->Write();
    h_Muon_Muon_invariant_mass->Write();
    h_Electron_Electron_invariant_mass->Write();
    h_Muon_Electron_invariant_mass_weighted->Write();
    h_Muon_Muon_invariant_mass_weighted->Write();
    h_Electron_Electron_invariant_mass_weighted->Write();
    h_leading_lepton_pt->Write();
    h_leading_lepton_pt_weighted->Write();
    h_LooseJets->Write();
    h_MediumJets->Write();
    h_TightJets->Write();
	h_dR_allJets->Write();
	h_dR_lbJets->Write();
	h_dR_mbJets ->Write();
      h_Apl_allJets->Write();
	h_Apl_lbJets->Write();
	h_Apl_mbJets->Write();
	h_Phi_allJets->Write();
	h_Phi_lbJets->Write();
	h_Phi_mbJets->Write();
	h_acopla_emu->Write();
	h_NJets->Write();
  h_vsPTandEta_onlye->Write();
  h_vsPTandEta_onlymu->Write();
  h_vsPTandEta_muande->Write();

h_Trigger->Write();
h_mu_3dsig->Write();
h_mu_3d ->Write();
h_mu_dxy->Write();
h_e_3dsig->Write();
h_e_3d->Write();
h_e_dxy->Write();

b_pt->Write();
jethole->Write();
ehole->Write();
}
#endif // Auxiliary_cpp
