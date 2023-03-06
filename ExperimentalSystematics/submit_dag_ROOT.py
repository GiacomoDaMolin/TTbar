from argparse import ArgumentParser
from subprocess import Popen, PIPE
import sys
import os
import json


def make_parser():
    parser = ArgumentParser(description="Submit Dag Jobs automatically.")
    parser.add_argument('-i', '--input_dir', type=str,
                        help="/afs/ base dir to submit from", default=None)
    parser.add_argument('-o', "--output_dir", type=str,
                        help="/eos/ base dir to store data", default=None)
    parser.add_argument('-e', "--executable", type=str, default=None,
                        help="full path to executable to run for \
each dag entry.")
    parser.add_argument('-j', "--json", type=str,
                        help="json file with datasets and xs..", default=None)
    parser.add_argument('-pr', "--proxy", type=str,
                        help="full path to voms cms user proxy", default=None)
    parser.add_argument('-s', '--systematics', type=bool,
                        help="If true, outputs also syst variations")
    return parser


def return_subfile(input_dir, base_dir, executable):
    if not executable.startswith("/"):
        executable = '/afs/cern.ch/user/g/gdamolin/Johan/\
TTbar/ExperimentalSystematics/'+executable
    arguments = 'Arguments = -i $(INFILE) -o $(OUTFILE) \
-x $(XS) -l $(LUMI) -d $(DATA) -s $(SYSTS) -p $(PROCESS)'

    file_str = f"basedir={input_dir}\n\
\n\
executable={executable}\n\
should_transfer_files = YES\n\
when_to_transfer_output = ON_EXIT\n\
\n\
output                = {base_dir}/out/$(ClusterId).$(ProcId).out\n\
error                 = {base_dir}/err/$(ClusterId).$(ProcId).err\n\
log                   = {base_dir}/log/$(ClusterId).$(ProcId).log\n\
\n\
+JobFlavour = \"workday\"\n\
{arguments}\n\
queue"
    return file_str


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


def write_dag(dagfile, subfile: str,
              infile: str, outfile: str, process:str,
              xs: float = None,
              lumi: int = None,
              data: bool = False,
              systs: bool= False,
              ):
    jobid = infile.split('/')[-1]
    infile = f"root://cms-xrd-global.cern.ch//{infile}"
    print(f"JOB {jobid} {subfile}", file=dagfile)
    print(f"VARS {jobid} INFILE=\"{infile}\" \
OUTFILE=\"{outfile}\" XS=\"{xs}\" LUMI=\"{lumi}\" \
DATA=\"{data}\" SYSTS=\"{systs} \" PROCESS=\"{process}\" ", file=dagfile)


def main():
    parser = make_parser()
    args = parser.parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir
    executable = args.executable
    systs=args.systematics
    for fd in [input_dir, output_dir, executable]:
        if not os.path.exists(fd):
            raise ValueError(f"{fd} does not exist.\
    Check your input or create it first.")
    data_json = open(args.json)
    dataS = json.load(data_json)
    for sample in dataS:
        basedir = input_dir+f"/{sample}"
        submit_file_str = return_subfile(input_dir=input_dir,
                                         base_dir=basedir,
                                         executable=executable,)
        for io_dir in [basedir, f"{output_dir}/{sample}"]:
            if not os.path.exists(io_dir):
                print(f"{io_dir} does not exist. Creating it now")
                os.mkdir(io_dir)
        dataset = dataS[sample]['dataset']
        
        xsec = dataS[sample]['xs']
        lumi = 59.82
        data = dataS[sample]['data']
        process=dataS[sample]['process']
        datafiles = run_dasgoclient(dataset=dataset)
        with open(f"{basedir}/{sample}.dag", 'w') as dagfile:
        	for file in datafiles:
        		write_dag(dagfile=dagfile,
                              subfile=f"{basedir}/{sample}.submit",
                              infile=file, outfile=f"{output_dir}/{sample}",
                              data=data, systs=systs, process=process,
                              xs=xsec, lumi=lumi)

                 
        with open(f"{basedir}/{sample}.submit", 'w') as file:
            print(submit_file_str, file=file)

        for directory in ['err', 'log', 'out']:
            if not os.path.exists(f"{basedir}/{directory}"):
                os.mkdir(f"{basedir}/{directory}")


if __name__ == "__main__":
    main()
