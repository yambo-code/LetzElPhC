#!/usr/bin/env python3
from driver import test_driver
import argparse
from test_list import tests
import os

### Test suite runs from here 
parser = argparse.ArgumentParser(description='LetzElPhC testsuite driver.')
parser.add_argument('-l','--lelphc', help='path to lelphc executable. Default is ./lelphc', required=False, metavar='', default='./lelphc')
parser.add_argument('-m','--mpirun', help='mpirun command. Default is mpirun', required=False, metavar='', default='mpirun')
parser.add_argument('-n','--ncpus', help='number of CPUs to use. Default is 1', required=False, metavar='', default=1, type=int)
parser.add_argument('-d','--download', help='Download test data.', action="store_true")
args = vars(parser.parse_args())


## if download requested, download the SAVE data
if args['download']:
    cwd = os.getcwd()
    os.chdir(os.path.join(cwd, 'dft_qe'))
    os.system("git init")
    os.system("git remote add origin https://github.com/muralidhar-nalabothula/Lelphc-Test-Data")
    os.system("git pull")
    os.system("git checkout main -f")
    os.system("git branch --set-upstream-to origin/main")
    os.chdir(cwd)
else :
    test_driver(tests, lelphc_cmd=args['lelphc'], mpirun_cmd=args['mpirun'], ncpus = args['ncpus'])