import subprocess
import uproot
from argparse import ArgumentParser
from subprocess import Popen, PIPE


def make_parser():
    parser = ArgumentParser(
        description="hadd files that have at least one entry.")
    parser.add_argument('-o', '--output_file', type=str,
                        help="name of output .root file", )
    parser.add_argument('-f', '--files', type=str,
                        nargs="*", help="input .root files")
    return parser


def run_hadd(outfile: str, infiles: str):
    """
    runs the hadd command
    """
    command = f"hadd {outfile} {infiles}"
    print(f"running:\n{command}\n")
    prcs = subprocess.Popen(command, shell=True,
                            stdout=PIPE, stderr=PIPE, encoding='utf-8')
    out, err = prcs.communicate()
    if err:
        print(err)
    else:
        print(out)


def get_branch


def main():
    parser = make_parser()
    args = parser.parse_args()
    filestring = ''
    for file in args.files:
        f = uproot.open(file)
        tree_keys = f['tout'].keys()
        m_pt_len = len(f['tout']['Muon_pt'].array())
        for key in tree_keys:
            keylen = len(f['tout'][key].array())
            if keylen != m_pt_len:
                raise ValueError(f"length of {key} does not coincide with len of Muon_pt")
        f.close()
        if m_pt_len > 0:
            filestring += file+" "
    run_hadd(outfile=args.output_file, infiles=filestring.strip())


if __name__ == "__main__":
    main()
            
            
