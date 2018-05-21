# This script defines a test case which computes one or more physical
# properties with a given model
#
# INPUTS:
#   model.calculator -- an ase.calculator.Calculator instance
#     this script can assume the calculator is checkpointed.
#
# OUTPUTS:
#   properties -- dictionary of key/value pairs corresponding
#     to physical quantities computed by this test

# standard ASE structure generation routines
import ase.io, os, sys

# set of utility routines specific this this model/testing framework
# the current model
import model, utilities

force_component_errors = []
ats = ase.io.read(os.path.join(os.path.dirname(__file__),'gp_iter6_sparse9k.xml.xyz'), index=':', format='extxyz')
for at in ats:
    if len(at) > 1:
        at.wrap()
        # ase.io.write(sys.stdout, at, format="extxyz")
        # sys.stdout.flush()
        try:
            at = utilities.evaluate(at) 
            try:
                dft_f = at.arrays['dft_force']
            except:
                try:
                    dft_f = at.arrays['DFT_force']
                except:
                    pass
            f = at.get_forces()
            for i in range(len(at)):
                force_component_errors.append((dft_f[i,0], f[i,0]-dft_f[i,0]))
                force_component_errors.append((dft_f[i,1], f[i,1]-dft_f[i,1]))
                force_component_errors.append((dft_f[i,2], f[i,2]-dft_f[i,2]))
        except:
            pass
        print len(force_component_errors)
        sys.stdout.flush()

        if 'GAP_TESTING_MODEL_RELOAD' in os.environ:
            del model.calculator
            reload(model)

# dictionary of computed properties - this is output of this test, to
#   be compared with other models
properties = {'force_component_errors' : force_component_errors }
