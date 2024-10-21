#!/usr/bin/env python3
from runner import test_driver
import argparse

## Add your tests here [foldername, .ini_filename]
tests = (
[
    ['.', 'test.ini'],
])


### Test suite runs from here 
parser = argparse.ArgumentParser(description='LetzElPhC testsuite driver.')
parser.add_argument('-l','--lelphc', help='path to lelphc executable. Default is ./lelphc', required=False, metavar='', default='./lelphc')
parser.add_argument('-m','--mpirun', help='mpirun command. Default is mpirun', required=False, metavar='', default='mpirun')
parser.add_argument('-n','--ncpus', help='number of CPUs to use. Default is 1', required=False, metavar='', default=1, type=int)
args = vars(parser.parse_args())

## defaults

test_driver(tests, lelphc_cmd=args['lelphc'], mpirun_cmd=args['mpirun'], ncpus = args['ncpus'])