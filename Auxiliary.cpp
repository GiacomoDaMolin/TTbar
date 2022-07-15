#ifndef Auxiliary_cpp
#define Auxiliary_cpp
#include <iostream>
// inputs are: size of arrays GenID and GenPar
// the pointer at the start of the array of GenParticles_PDGID
// the pointer at the start of the array of GenParticles_Parent
// initialID: last index of the particle (ex: Electron_genPartIdx[j]=initialID)

// should be called from func as isFromW(MAX_ARRAY_SIZE,GenPart_pdgId,GenPart_genPartIdxMother,Electron_genPartIdx[j])

bool isFromW(int size, Int_t *GenId, Int_t *GenParent, int initialID)
{
	// retrieve first PDG ID number
	int startPdg = GenId[initialID];
	int newID = initialID, newPdg = startPdg;
	// look for the parent; if the parent is of same PDGID of starting particle, iterate until parent is different particle
	while (newPdg == startPdg)
	{
		if (newID > size)
		{
			std::cout << "WARNING: index " << newID << " exceeding max size " << size << endl;
		}
		newID = GenParent[newID];
		newPdg = GenId[newID];
		if (abs(newPdg) == 24)
			return true;
	}
	return false;
}


void printMCTree(int size, Int_t *GenId, Int_t *GenParent, Int_t initialID)
{
	// retrieve first PDG ID number
	Int_t startPdg = GenId[initialID]; 
	Int_t newPdg = startPdg;
	Int_t newID = initialID; 
	// look for the parent; if the parent is of same PDGID of starting particle, iterate until parent is different particle
	while (newPdg == startPdg)
	{
		if (newID > size)
		{
			std::cout << "WARNING: index " << newID << " exceeding max size " << size << endl;
		}
		newID = GenParent[newID];
		newPdg = GenId[newID];
		std::cout << "ID: " << newID << " PDG: " << newPdg << endl;
	}
}

#endif // Auxiliary_cpp