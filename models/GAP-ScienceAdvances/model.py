import os

from quippy import Potential
import __builtin__

orig_dir = os.getcwd()
model_dir = os.path.dirname(__file__)
if model_dir != '':
    os.chdir(model_dir)

if os.path.exists('gp_iter5_with111adatom_with3x3das_sparse2x.xml.sparseX.GAP_2017_5_20_60_4_23_20_5121.bz2'):
    os.system('bunzip2 gp_iter5_with111adatom_with3x3das_sparse2x.xml.sparseX.GAP_2017_5_20_60_4_23_20_5121.bz2')

try:
    if hasattr(__builtin__, 'mpi_glob'):
        calculator = Potential(init_args='Potential xml_label="GAP_2017_5_20_60_4_23_20_512"',
                                               param_filename='gp_iter5_with111adatom_with3x3das_sparse2x.xml', mpi_obj=mpi_glob)
    else:
        calculator = Potential(init_args='Potential xml_label="GAP_2017_5_20_60_4_23_20_512"',
                                               param_filename='gp_iter5_with111adatom_with3x3das_sparse2x.xml')
    Potential.__str__ = lambda self: '<GAP Potential>'
finally:
    os.chdir(orig_dir)

no_checkpoint = True

name = 'GAP-ScienceAdvances'
