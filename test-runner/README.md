# script to run a test using a model

Usage: `Usage: run-model-test.py [-h] [--no_redirect_stdout] [--force] model_name test_name`

Arguments:
* `--help` or `-h`: optional flag to disable writing output to `model-MODEL-test-TEST.txt`
* `--no-redirect-stdout` or `-n`: optional flag to disable writing output to `model-MODEL-test-TEST.txt`
* `--force` or `-f`: optional flag to run test even if it appears to have already been run
* `model_name`: name of model (path under `models/` directory)
* `test_name`: name of test (path under `tests/` directory)

Output consists of stdout contents in `model-MODEL-test-TEST.txt` (unless redirection is disabled), main test 
results in `model-MODEL-test-TEST-properties.json`, and a directory (`model-MODEL-test-TEST/`) which contains intermediate files, 
relaxed geometries, etc.

If you are using the `libatomsquip/quip` to run the tests, you need to `export LAMMPS_LIB=` before invoking `run-model-test.py`.