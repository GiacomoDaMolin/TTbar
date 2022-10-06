#!/afs/cern.ch/user/j/jowulff/miniconda3/envs/shep/bin/python
import uproot
import numpy as np
from argparse import ArgumentParser


def make_parser():
    parser = ArgumentParser(description="Divide an hadded sample by the sum of\
the sum of weights of each file.")
    parser.add_argument('-i', '--input_files', type=str, default=None,
                        help="Input .root file(s)", nargs='*')
    return parser


def normalise(filename: str,):
    file = uproot.open(filename)
    sum_w = np.sum(file['Run_out']['SumW'].array())
    weights = file['tout']['weight'].array()
    gen_weights = file['tout']['genWeight'].array()
    weights = weights / sum_w
    gen_weights = gen_weights/sum_w
    outfile = uproot.update(filename)
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
