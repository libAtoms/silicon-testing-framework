# parameters from AvD e-mail 1/18/2017, as used in
# Buehler, M.J., van Duin, A.C.T. and Goddard, W.A. (2006) Multiparadigm modeling of dynamical crack propagation in silicon using a reactive force field. Physical Review Letters 96, 095505.

import os
import sys

try:
    from ase.calculators.lammpslib import LAMMPSlib
except:
    sys.exit("Failed to import lammpslib module")

model_dir = os.path.dirname(__file__)

header=['units real',
        'atom_style charge',
        'atom_modify map array sort 0 0']

cmds = ["pair_style reax/c NULL checkqeq no",
        "pair_coeff * * {}/ffield_APL_2014_SiO Si".format(model_dir)]

try:
   calculator = LAMMPSlib(lmpcmds = cmds,keep_alive=True,lammps_header=header,atom_types={"Si":1}, lammps_name=os.environ['LAMMPS_NAME'])
except:
   raise RuntimeError

no_checkpoint = True

name = 'ReaxFF'
