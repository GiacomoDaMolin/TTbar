from argparse import ArgumentParser
from asyncio import subprocess
from glob import glob
from subprocess import Popen


def make_parser():
    parser = ArgumentParser(description="Resubmit dags that created\
a rescue file. Run this script on a directory and it will check it's\
subdirectories for rescue00N files where N is an int.")
    parser.add_argument('-d', '--directory', type=str, help='parent\
directory to check subdirectories from.')
    parser.add_argument('-n', '--number', type=int, help='the number added\
at the end of the rescue dag file.')
    return parser


def resubmit(directory: str, n: int):
    if directory.endswith('/'):
        directory = directory[:-1]
    rescuefiles = glob(directory+f"/*/*.rescue00{n}")
    print(f"Found {len(rescuefiles)} .rescue00{n}")
    for rfile in rescuefiles:
        dagfile = '.'.join(rfile.split('.')[:-1])
        cmd = f"condor_submit_dag {dagfile}"
        print(f"running: {cmd}")
        proc = Popen(cmd, shell=True,
                     stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = proc.communicate()
        if err:
            print(err)
        else:
            print("success")


def main():
    parser = make_parser()
    args = parser.parse_args()
    resubmit(directory=args.directory, n=args.number)


if __name__ == "__main__":
    main()
