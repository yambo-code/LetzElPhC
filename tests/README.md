## How to run the test suite?
The test suite is written in ```python3```. To run the test suite, the following modules are required: <br/>
1) Numpy <br/>
2) Scipy <br/>
3) netCDF4 <br/>

```console
python3 run_tests.py --ncpus=1 --mpirun=mpirun --lelphc=./lelphc
## For more information about the available options use the below command
python3 run_tests.py --help 
``` 

## Adding new tests to the test suite
Adding tests is straightforward. Append a two-element list to the ```tests``` list located in ```test_list.py```. 
The list should follow this format: ```[test_folder, test_ini_file]```, where ```test_folder``` is the relative 
path to the test directory. ```test_ini_file``` is the ```.ini``` configuration file that must be present in the 
```test_folder```. The ```test_ini_file``` defines the list of tests in the test folder. It follows this structure:


```dosini
; Tests are defined using sections. Each section defines a test.
; The DEFAULT section (all uppercase) contains inputs that apply to all other sections.
;
; The ndb.elph and ndb.Dmat reference files are set with elph_db_ref and dmat_db_ref, respectively.
; Section names follow the structure: [testname code], where code is typically qe. 
; If only one name is provided (i.e., [testname]), the code is assumed to be qe.
;
[DEFAULT] ; 
start_bnd       = 1 
end_bnd         = 40 
save_dir        = /Users/murali/phd/one_phonon_raman/wse2/SAVE  
ph_save_dir     = /Users/murali/phd/one_phonon_raman/wse2/ph_save  
kernel          = dfpt
elph_db_ref     =  ../ndb.elph_ref_test
dmat_db_ref     =  ../ndb.Dmats_ref_test
;
[test1 qe]
convention      = yambo
;
[test2]
convention      = standard
```



