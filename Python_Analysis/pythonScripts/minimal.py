#!/afs/cern.ch/user/j/jowulff/miniconda3/envs/shep/bin/python
import uproot
import numpy as np
import awkward as ak

from argparse import ArgumentParser


def make_parser():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", type=str, 
                        default=None, help='Input Root File')
    parser.add_argument("-o", "--output", type=str,
                        default=None, help='Output File')
    return parser


def skimming(filename, ofilename, xs=None, lumi=None, mc_flag=False,):
    outfile = uproot.recreate(ofilename)
    out_dict = {'Muon_pt': 'var *float64'}
    outfile.mktree("tout", out_dict)
    filter_names = ['Muon_pt', ]
    #file_handler = uproot.MultithreadedFileSource
    #file = uproot.open(filename, file_handler=file_handler)
    file = uproot.open(filename,)
    tree = file['Events']
    trigger_cut = "HLT_IsoMu24 | HLT_Ele32_WPTight_Gsf"
    for events in tree.iterate(
            filter_name=filter_names, cut=trigger_cut,
            entry_stop=10000):
        tout_dict = {'Muon_pt': events['Muon_pt'],}
        outfile['tout'].extend(tout_dict)


def main():
    parser = make_parser()
    args = parser.parse_args()
    skimming(filename=args.input,
             ofilename=args.output,)


if __name__ == "__main__":
    main()
