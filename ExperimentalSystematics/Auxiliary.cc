#ifndef Auxiliary_cpp
#define Auxiliary_cpp
#include <iostream>
#include<vector>
#include<string>
using std::vector;

float getTopPtWeight(Int_t * pdgId,Int_t *statusFlags,Float_t * pt, Int_t Ngen) {
    float wgt=1.0;
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

TH1D* cloneDims1d(TH1D* hist, string newname){
  if (hist == NULL) return NULL;
  TH1D* cloneHist = new TH1D(Form("%s_%s", hist->GetName(),newname), 
                             Form("%s_%s", hist->GetName(),newname), 
                             hist->GetNbinsX(),
                             hist->GetXaxis()->GetXmin(),hist->GetXaxis()->GetXmax());
                             
  cloneHist->Sumw2();
  return cloneHist;
}

/*size_t found = ofile.find_last_of("/");
    string oname=ofile.substr(found+1);
    string path=ofile.substr(0,found);
    string Tauname=path+"/Tau_"+oname;*/

void HistWrite(){
}
#endif // Auxiliary_cpp
