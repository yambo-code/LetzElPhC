## How to run test suite ?
The test suite is written in python3. To run the test suite, we need the following modules <br/>
1) Numpy <br/>
2) Scipy <br/>
3) netCDF4 <br/>

```console
python3 tests.py --ncpus=1 --mpirun=mpirun --lelphc=./lelphc
```
## Adding new tests to test suite
Adding tests is straight forward. Append a two element list to ```tests``` list present in ```tests.py```. The list should have the following format: [foldername, test_ini_file]



