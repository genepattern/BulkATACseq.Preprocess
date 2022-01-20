# Input BAM File, optional control file

#!/usr/bin/env python3

import sys
import argparse
from io import StringIO
import os
import subprocess
from subprocess import PIPE, STDOUT


### Code for Inputs
parser = argparse.ArgumentParser()

### Module Required Arguments
parser.add_argument("-t", "--treatment",
                    type = str,
                    help ="Name of the file to be read")

parser.add_argument('-k',"--cutoff", help='Alignment number cutoff') 
parser.add_argument('-p','--paired', help='Data is paired-end')


args = parser.parse_args()

buff = StringIO()

buff.write("assign_multimappers.py -t")

file_list = args.treatment

buff.write(" ")
buff.write(file_list)

if args.cutoff:
    buff.write(" -k ")
    buff.write(args.cutoff)

if args.paired:
    buff.write(" --paired-end")

command_str = buff.getvalue()
print(command_str)

# subprocess = subprocess.Popen(command_str, shell = True, stdout=PIPE)
# subprocess_return = subprocess.stdout.read()
# print(subprocess_return)
