#include<iostream>
#include<string>

void upd(string namefile) { 
TFile *f = new TFile(namefile.c_str(),"update"); 
TTree *T = (TTree*)f->Get("tout");
TTree *W = (TTree*)f->Get("weights");
double nW; 

W->SetBranchStatus("*",0);
W->SetBranchStatus("normalised_weights",1);
W->SetBranchAddress("normalised_weights",&nW); 
TBranch *bout = T->Branch("norm_W",&nW,"norm_W/D"); 
Long64_t nentries = T->GetEntries(); 
Long64_t n_entries_check= W->GetEntries(); 
if (nentries!=n_entries_check) {
	std::cout<<"Different number of entries in trees"<<std::endl;
	return;
	}
 for (Long64_t i=0;i<nentries;i++) { 
 	W->GetEntry(i);
 	T->GetEntry(i);
 	bout->Fill(); 
 	} 
 	
T->Print(); 
T->Write(); 
delete f; 
}
