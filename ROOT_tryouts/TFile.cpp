#include <iostream>

#include "TFile.h"
#include "TTree.h"

void TFileExample(string inputfile, string outputfile)
{
    TFile *file = new TFile(inputfile.c_str());
    // create the output file
    TFile *outfile = new TFile(outputfile.c_str(), "RECREATE");
    TTree *tree = (TTree *)file->Get("tree");
    tree->Print();
    file->Close();
}

int main(int argc, char **argv)
{
    string inputFile = argv[1];
    string outputFile = argv[2];
    TFileExample(inputFile, outputFile);
}
