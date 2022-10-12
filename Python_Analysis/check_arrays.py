import uproot
from argparse import ArgumentParser


def make_parser():
    parser = ArgumentParser(description="Check files for good TBranches.")
    parser.add_argument('-f', '--files', type=str,
                        nargs="*", help="input .root files")
    return parser


def main():
    parser = make_parser()
    args = parser.parse_args()
    for file in args.files:
        f = uproot.open(file)
        tree_keys = f['tout'].keys()
        arr_lengths = set([]) 
        for idx, key in enumerate(tree_keys):
            try:
                arr = f['tout'][key].array()
            except:
                print(f"{key} not accesible for {file}")

            if len(arr) not in arr_lengths:
                if idx > 0:
                    print(f"Length of {key} is: {len(arr)}")
                    print(f"Length of {tree_keys[idx-1]} is: \
{len(f['tout'][tree_keys[idx-1]].array())}")
                arr_lengths.add(len(arr))
        f.close()
        if len(arr_lengths) == 1:
            print(f"All branches of equal length for {file}")


if __name__ == "__main__":
    main()
