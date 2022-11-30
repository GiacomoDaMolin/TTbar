#!/eos/user/j/jowulff/miniconda3/bin/python
import json
import uproot
import numpy as np
import sys
from argparse import ArgumentParser
from subprocess import Popen, PIPE


def make_parser():
    parser = ArgumentParser(
        description='Get the Sum of the Weights for a MC dataset')
    parser.add_argument('-j', '--json_file', type=str, help="input json file")
    return parser


def run_dasgoclient(dataset: str):
    """
    Runs dasgoclient and returns a list of files for a given dataset
    """
    cmd = f'dasgoclient -query="file dataset={dataset}"'
    process = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, encoding='utf-8')
    out, err = process.communicate()
    if err:
        print(err)
        sys.exit(1)
    else:
        return out.split()


def main():
    parser = make_parser()
    args = parser.parse_args()
    with open(args.json_file, 'r+') as json_file:
        file_data = json.load(json_file)
        for sample in file_data:
            print(sample)
            if 'Sum_w' in file_data[sample].keys():
                if file_data[sample]['Sum_w'] > 1.:
                    print(f"{sample} has key Sum_w with value {file_data[sample]['Sum_w']}")
                    continue
            print("hello")
            print(f"Getting weights for {file_data[sample]['dataset']}")
            files = run_dasgoclient(dataset=file_data[sample]['dataset'])
            ne, sum_w = [], []
            xrootd_string = 'root://cms-xrd-global.cern.ch//'
            for n, file in enumerate(files):
                print(
                    f"opening file number {n+1} out of {len(files)}", end='\r')
                while True:
                    try:
                        f = uproot.open(xrootd_string+file)
                        ne.append(f['Runs']['genEventCount'].array())
                        sum_w.append(f['Runs']['genEventSumw'].array())
                        f.close()
                    except OSError:
                        continue
                    break

            json_dict = {"Sum_w": float(
                np.sum(sum_w)), "N_events": int(np.sum(ne))}
            print(json_dict)
            file_data[sample].update(json_dict)
            json_file.seek(0)
            json.dump(file_data, json_file, indent=4)


if __name__ == '__main__':
    main()
