#!/usr/bin/env python


import os
import sys
import datetime
import json
import time

import __builtin__

if 'USE_MPI' in os.environ:
    import mpi4py
    from quippy.mpi_context import MPI_context
    __builtin__.mpi_glob = MPI_context()
    __builtin__.do_io = (mpi4py.MPI.COMM_WORLD.Get_rank() == 0)
    print "rank ",mpi4py.MPI.COMM_WORLD.Get_rank(), "io", __builtin__.do_io
else:
    __builtin__.do_io = True

try:
    from ase.calculators.checkpoint import CheckpointCalculator
except:
    pass

import argparse

import argparse
parser = argparse.ArgumentParser(description='run a test with a model')
parser.add_argument('--no_redirect_stdout','-n',dest='redirect_stdout',action='store_false', help='do not redirect output to file')
parser.add_argument('--force','-f',action='store_true', help='run even if already seems to be complete (*-properties.json file exists')
parser.add_argument('model_name',action='store', type=str, help='model directory')
parser.add_argument('test_name',action='store', type=str, help='test directory')
args = parser.parse_args()

dir_name = 'model-{0}-test-{1}'.format(args.model_name, args.test_name)
if not os.path.exists(dir_name) and __builtin__.do_io:
    os.mkdir(dir_name)
time.sleep(1)
os.chdir(dir_name)

json_file_name = os.path.join('..', 'model-{0}-test-{1}-properties.json'.format(args.model_name, args.test_name))
if not args.force and os.path.isfile(json_file_name) and os.path.getsize(json_file_name) > 0:
    print "%s already exists and is not empty, not rerunning test" % json_file_name
    sys.exit(0)

share_dir = os.path.join('..', '..', 'share')
model_dir = os.path.join('..', '..', 'models', args.model_name)
test_dir = os.path.join('..', '..', 'tests', args.test_name)

sys.path.insert(0, share_dir)
sys.path.insert(0, model_dir)
sys.path.insert(0, test_dir)

if args.redirect_stdout:
    _stdout, _stderr = sys.stdout, sys.stderr
    if __builtin__.do_io:
        log = open('../model-{0}-test-{1}.txt'.format(args.model_name, args.test_name), 'w', 0)
        sys.stdout, sys.stderr = log, log
    else:
        sys.stdout = open(os.devnull, "w")
        sys.stderr = open(os.devnull, "w")

import logging
logging.basicConfig(stream=sys.stdout, level=logging.INFO)
logger = logging.getLogger('ase.optimize.precon')
logger.propagate = True
logger.setLevel(logging.INFO)

print 'Model {0}, Test {1}'.format(args.model_name, args.test_name)
print 'Test run at {0}'.format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M"))
print

model_file = os.path.join(model_dir, 'model.py')
print 'model file:',model_file
print '='*60
sys.stdout.write(open(model_file).read())
print '='*60

test_file = os.path.join(test_dir, 'test.py')
print 'test file:', test_file
print '='*60
sys.stdout.write(open(test_file).read())
print '='*60

import model # import and run the current model

# adapt model to test, e.g. to set output directory name
if hasattr(model, 'start'):
    model.start(args.test_name)

# create a checkpoint file for this specific (model, test) combination
if ('SI_GAP_TEST_CHECKPOINT' in os.environ and os.environ['SI_GAP_TEST_CHECKPOINT'] != "") or not hasattr(model, 'no_checkpoint') or not model.no_checkpoint:
   checkpoint_file = 'model-{0}-test-{1}.db'.format(args.model_name, args.test_name)
   model.calculator = CheckpointCalculator(model.calculator, db=checkpoint_file)

import test  # import and run the current test

print '='*60
print 'Property calculation output:'
print

# serialise results in machine readable JSON format
json_file = open(json_file_name, 'write')
json.dump(test.properties, json_file)
json_file.close()

print
print 'Summary of computed properties:'
print test.properties

print '='*60

if args.redirect_stdout:
    sys.stdout, sys.stderr = _stdout, _stderr
    if __builtin__.do_io:
        log.close()

if hasattr(model, 'shutdown'):
    model.shutdown()
