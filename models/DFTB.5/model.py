import os

from quippy import Potential
import __builtin__

orig_dir = os.getcwd()
model_dir = os.path.dirname(__file__)
os.chdir(model_dir)
try:
    if hasattr(__builtin__, 'mpi_glob'):
        calculator = Potential(init_args='TB DFTB use_k_density k_density=20.0 Fermi_T=0.05',
                                               param_filename='tightbind.parms.DFTB.pbc-0-1.xml', mpi_obj=mpi_glob)
    else:
        calculator = Potential(init_args='TB DFTB use_k_density k_density=20.0 Fermi_T=0.05',
                                               param_filename='tightbind.parms.DFTB.pbc-0-1.xml')
finally:
    os.chdir(orig_dir)

no_checkpoint = True

name = 'DFTB.5'
