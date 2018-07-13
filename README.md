# silicon-testing-framework
Testing framework for silicon interatomic potentials

This is a lightweight set of python scripts that facilitate the testing of predictions of interatomic potentials for a wide range of material properties. The [Atomic Simulation Environment](https://wiki.fysik.dtu.dk/ase) is used to "tie together" the potentials and the tests. ASE includes an interface to [LAMMPS](http://lammps.sandia.gov), so models can be created through that route as well. 

Feel free to fork this code to create tests for other materials! 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1250555.svg)](https://doi.org/10.5281/zenodo.1250555)

* `models` directory containing each model (DFTB and analytical interatomic potentials)
* `tests` scripts for doing tests
* `share` utilities used by tests
* `test-runner` place to run tests, including script "run_model_test.py" that actually runs a given test with a given model 

