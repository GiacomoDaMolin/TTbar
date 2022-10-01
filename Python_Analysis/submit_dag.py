from argparse import ArgumentParser
from asyncio import subprocess
from subprocess import Popen
import os
import json


def make_parser():
    parser = ArgumentParser(description="Submit Dag Jobs automatically.")
    parser.add_argument('-i', '--input_dir', type=str,
                        help="/afs/ base dir to submit from", default=None)
    parser.add_argument('-o', "--output_dir", type=str,
                        help="/eos/ base dir to store data", default=None)
    parser.add_argument('-e', "--executable", type=str, default=None,
                        help="executable to run for each dag entry.")
    parser.add_argument('-j', "--json", type=str,
                        help="json file with datasets and xs..", default=None)
    parser.add_argument('-p', "--proxy", type=str,
                        help="voms cms user proxy", default=None)
    return parser


def return_subfile(input_dir, base_dir, executable):
    file_str = f"basedir={input_dir}\n\
\n\
executable={input_dir}/{executable}\n\
should_transfer_files = YES\n\
when_to_transfer_output = ON_EXIT\n\
\n\
output                = {base_dir}/out/hello.$(ClusterId).$(ProcId).out\n\
error                 = {base_dir}/err/hello.$(ClusterId).$(ProcId).err\n\
log                   = {base_dir}/log/hello.$(ClusterId).log\n\
\n\
+JobFlavour = \"longlunch\"\n\
Arguments = -i $(INFILE) -o $(OUTFILE) -x $(XS) -l $(LUMI) -p $(PROXY)\n\
queue"
    return file_str


def main():
    parser = make_parser()
    args = parser.parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    executable = args.executable
    proxy = args.proxy
    for fd in [input_dir, output_dir, input_dir+'/'+proxy,
               input_dir+'/'+executable]:
        if not os.path.exists(fd):
            raise ValueError(f"{fd} does not exist.\
    Check your input or create it first.")
    data_json = open(args.json)
    data = json.load(data_json)
    jobscript = '/eos/user/j/jowulff/TTbar/PythonAnalysis/write_dag.sh'
    for sample in data:
        basedir = input_dir+f"/{sample}"
        submit_file_str = return_subfile(input_dir=input_dir,
                                         base_dir=basedir,
                                         executable=executable)
        for io_dir in [basedir, f"{output_dir}/{sample}"]:
            if not os.path.exists(io_dir):
                print(f"{io_dir} does not exist. Creating it now")
                os.mkdir(io_dir)
        dataset = data[sample]['dataset']
        xsec = data[sample]['xs']
        lumi = 59.82
        cmd = f"{jobscript} -f {basedir}/{sample}.submit \
-j {basedir}/{sample}.dag -d {dataset} -o {output_dir}/{sample} \
-x {xsec} -l {lumi} -p {input_dir+'/'+proxy}"
        write_dagfile = Popen(cmd.split(), stdout=subprocess.PIPE)
        wd_out, wd_err = write_dagfile.communicate()
        with open(f"{basedir}/{sample}.submit", 'w') as file:
            print(submit_file_str, file=file)

        for directory in ['err', 'log', 'out']:
            if not os.path.exists(f"{basedir}/{directory}"):
                os.mkdir(f"{basedir}/{directory}")


if __name__ == "__main__":
    main()
