#!/afs/cern.ch/user/j/jowulff/miniconda3/envs/shep/bin/python
import uproot
import numpy as np
from argparse import ArgumentParser
import hist


def make_parser():
    parser = ArgumentParser(description="Divide an hadded sample by the sum of\
the sum of weights of each file.")
    parser.add_argument('-i', '--input_files', type=str, default=None,
                        help="Input .root file(s)", nargs='*')
    return parser


def normalise(filename: str,):
    h_Muon_pt_weighted = hist.Hist(hist.axis.Regular(
        bins=100, start=0, stop=200, name="Muon pt"))
    h_Muon_eta_weighted = hist.Hist(hist.axis.Regular(
        bins=100, start=-5, stop=5, name="Muon eta"))
    h_Electron_pt_weighted = hist.Hist(hist.axis.Regular(
        bins=100, start=0, stop=200, name="Electron pt"))
    h_Electron_eta_weighted = hist.Hist(hist.axis.Regular(
        bins=100, start=-5, stop=5, name="Electron eta"))
    h_Muon_Electron_invariant_mass_weighted = hist.Hist(hist.axis.Regular(
        bins=100, start=12, stop=412, name="Muon Electron Inv. Mass"))
    h_leading_lepton_pt_weighted = hist.Hist(hist.axis.Regular(
        bins=45, start=20, stop=200, name="leading lepton pt"))

    file = uproot.open(filename)
    sum_w = np.sum(file['Run_out']['SumW'].array())
    weights = file['tout']['weight'].array()
    gen_weights = file['tout']['genWeight'].array()
    weights = weights / sum_w
    gen_weights = gen_weights/sum_w
    mu_e_inv_mass = file['tout']['mu_e_inv_mass']
    leading_lepton_pt = file['tout']['leading_lepton_pt']
    muon_pt = file['tout']['muon_pt']
    muon_eta = file['tout']['muon_eta']
    electron_pt = file['tout']['electron_pt']
    electron_eta = file['tout']['electron_eta']
    file.close()
    h_Muon_Electron_invariant_mass_weighted.fill(
        mu_e_inv_mass, weight=weights)
    h_leading_lepton_pt_weighted.fill( leading_lepton_pt, weight=weights)
    h_Muon_eta_weighted.fill(muon_eta, weight=weights)
    h_Muon_pt_weighted.fill(muon_pt, weight=weights)
    h_Electron_pt_weighted.fill(electron_pt, weight=weights)
    h_Electron_eta_weighted.fill(electron_eta, weight=weights)
    outfile = uproot.update(filename)
    outfile['h_Muon_pt_weighted'] = h_Muon_pt_weighted
    outfile['h_Muon_eta_weighted'] = h_Muon_eta_weighted
    outfile['h_Electron_pt_weighted'] = h_Electron_pt_weighted
    outfile['h_Electron_eta_weighted'] = h_Electron_eta_weighted
    outfile['h_Muon_Electron_invariant_mass_weighted'] = h_Muon_Electron_invariant_mass_weighted
    outfile['h_leading_lepton_pt_weighted'] = h_leading_lepton_pt_weighted
    outfile['weights'] = {"normalised_weights": weights,
                          "normalised_gen_weights": gen_weights}


def main():
    parser = make_parser()
    args = parser.parse_args()
    infiles = args.input_files
    for file in infiles:
        normalise(file)


if __name__ == '__main__':
    main()
