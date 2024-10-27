# This file is part of test suite.
# This file contains the following functions 
# 1) Convert the non-portable binary file to portable numpy and vice versa
# 2) Generate float/double SAVE databases (cdo -b f64 copy input.nc output.nc)


import numpy as np
from netCDF4 import Dataset
import os
import glob


def binary2npy(bin_fine,read_format=np.double,write_format=np.single):
    ## Generates a portable .npy file from C/Fortran binary file.
    # read_format is datatype in binary format
    # write_format datatype used in dumped .npy files
    npy_file_name = 'npy_'+bin_fine.strip()+'.npy'
    tmp_data = np.ascontiguousarray(np.fromfile(bin_fine, dtype=read_format).astype(write_format))
    np.save(npy_file_name,tmp_data)


def npy2binary(bin_fine,write_format=np.double):
    ## Generates a C/Fortran binary file from portable .npy file
    ## write_format is datatype used in dumped binary file
    npy_file_name = 'npy_'+bin_fine.strip()+'.npy'
    tmp_data = np.ascontiguousarray(np.load(npy_file_name).astype(write_format))
    tmp_data.tofile(bin_fine)


def generate_portable_ph_save(ph_save_dir,remove=True):
    ## converts all dvscf* files in ph_save_dir to numpy format
    ## if remove = True, will remove non-portable dvscf files
    cwd = os.getcwd()
    os.chdir(os.path.join(cwd, ph_save_dir))
    dvscfs_names = glob.glob('dvscf*')

    for idvscf in dvscfs_names:
        binary2npy(idvscf)

    if remove:
        os.system('rm dvscf* > /dev/null 2>&1')
    
    os.chdir(cwd)


def generate_binary_ph_save(ph_save_dir,remove=False):
    ## converts all npy_dvscf* files in ph_save_dir to binary format
    ## if remove = True, will remove  numpy dvscfs files
    cwd = os.getcwd()
    os.chdir(os.path.join(cwd, ph_save_dir))
    dvscfs_names = glob.glob('npy_dvscf*')

    for idvscf in dvscfs_names:
        binary2npy(idvscf.strip().replace('npy_','').replace('.npy','').strip())

    if remove:
        os.system('rm npy_dvscf* > /dev/null 2>&1')
    
    os.chdir(cwd)


