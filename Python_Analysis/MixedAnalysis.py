import uproot
import numpy as np
import awkward as ak

import vector

import hist
from argparse import ArgumentParser


def make_parser():
    parser = ArgumentParser()
    parser.add_argument("-i", "--input", nargs='*',
                        default=None, help='Input Root File')
    parser.add_argument("-o", "--output", nargs=1,
                        default=None, help='Output File')
    parser.add_argument("-x", "--cross_section", nargs=1,
                        default=-1, help="Cross Section of MC Dataset")
    parser.add_argument("-l", "--int_luminosity", nargs=1,
                        default=-1, help="Integrated Lumi to scale to")
    parser.add_argument("-m", "--mc", nargs=1,
                        action='store_true', help='MC Flag. If set, mc==True')
    parser.add_argument("-s", "--signal", nargs=1, action='store_true',
                        help='Signal Flag. If set, signal==True. Only allowed if mc == True.')
    return parser


def spacial_invert(vec):
    return -vec.to_Vector3D().to_Vector4D()+vector.obj(px=0, py=0, pz=0, E=vec.E)


filename = "./input/MC/Signal/019426EE-3D50-1249-B266-F6DBA0AFE3B5.root"

file = uproot.open(filename)
