#!/usr/bin/env python3


# This file is part of test suite.
# This file contains the following functions 
# 1) Convert the non-portable binary file to portable numpy and vice versa
# 2) Generate float/double SAVE databases.


import numpy as np
from netCDF4 import Dataset
import os
import glob


def binary2npy(bin_file,read_format=np.double,write_format=np.single):
    ## Generates a portable .npy file from C/Fortran binary file.
    # read_format is datatype in binary format
    # write_format datatype used in dumped .npy files
    npy_file_name = 'npy_'+bin_file.strip()+'.npy'
    tmp_data = np.ascontiguousarray(np.fromfile(bin_file, dtype=read_format).astype(write_format))
    np.save(npy_file_name,tmp_data)


def npy2binary(bin_file,write_format=np.double):
    ## Generates a C/Fortran binary file from portable .npy file
    ## write_format is datatype used in dumped binary file
    npy_file_name = 'npy_'+bin_file.strip()+'.npy'
    tmp_data = np.ascontiguousarray(np.load(npy_file_name).astype(write_format))
    tmp_data.tofile(bin_file)


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
    npy_dvscfs = glob.glob('npy_dvscf*')
    bin_dvscfs = glob.glob('dvscf*')

    if (len(npy_dvscfs) != 0 and len(bin_dvscfs) == 0):
        for idvscf in npy_dvscfs:
            npy2binary(idvscf.strip().replace('npy_','').replace('.npy','').strip())

        if remove:
            os.system('rm npy_dvscf* > /dev/null 2>&1')
    os.chdir(cwd)



def nc_convert_types(ncfile, dtype_in, dtype_out,replace):
    ## changes all variables from float to double in a given netcdf file
    ## dtype_in is list
    ## Adapted from https://stackoverflow.com/a/49592545
    out_file_tmp = ncfile.strip()+'_tmp'
    with Dataset(ncfile) as src, Dataset(out_file_tmp, "w") as dst:
        dst.setncatts(src.__dict__)
        # copy dimensions
        for name, dimension in src.dimensions.items():
            dst.createDimension(
                name, (len(dimension) if not dimension.isunlimited() else None))
        # copy all file data except for the excluded
        for name, variable in src.variables.items():
            if variable.datatype in dtype_in:
                x = dst.createVariable(name, dtype_out, variable.dimensions)
            else :
                x = dst.createVariable(name, variable.datatype, variable.dimensions)
            dst[name][:] = src[name][:]
            # copy variable attributes all at once via dictionary
            dst[name].setncatts(src[name].__dict__)

    if replace:
        os.system("rm %s > /dev/null 2>&1" %(ncfile))
        os.system("mv %s %s > /dev/null 2>&1"%(out_file_tmp,ncfile))


def convert_save_dbs(save_dir, dtype_in, dtype_out,replace=True):
    cwd = os.getcwd()
    os.chdir(os.path.join(cwd, save_dir))
    save_dbs = glob.glob('ns.*') + glob.glob('ndb.*')

    is_save_in_dtype_out = False

    ## check if save is in given dtype
    ns_db_test = Dataset('ns.db1','r')

    if ns_db_test['EIGENVALUES'].datatype == dtype_out:
        is_save_in_dtype_out = True
    
    ns_db_test.close()

    if not is_save_in_dtype_out:
        for idb in save_dbs:
            nc_convert_types(idb, dtype_in, dtype_out, replace)
    
    os.chdir(cwd)


## Use this as a program
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser()

    group = parser.add_mutually_exclusive_group()

    group.add_argument('-n','--to_npy', \
        help='ph_save directory. Converts binary dvscfs to portable npy databases.', \
        required=False, metavar='', default='')
    
    group.add_argument('-b','--to_binary', \
        help='ph_save directory. Converts portable npy databases to binary dvscf files.', \
            required=False, metavar='', default='')

    group.add_argument('-f','--to_float', \
        help='SAVE directory. Creates a SAVE as if it was created by single precision p2y.', \
            required=False, metavar='', default='')
    
    group.add_argument('-d','--to_double', \
        help='SAVE directory. Creates a SAVE as if it was created by double precision p2y.', \
            required=False, metavar='', default='')

    # parse all variables
    args = vars(parser.parse_args())

    if args['to_binary'] != '':
        generate_binary_ph_save(args['to_binary'].strip(),remove=False)
    elif args['to_npy'] != '':
        generate_portable_ph_save(args['to_npy'].strip(),remove=True)
    elif args['to_double'] != '':
        convert_save_dbs(args['to_double'].strip(), [np.double, np.single], np.double,replace=True)
    elif args['to_float'] != '':
        convert_save_dbs(args['to_float'].strip(), [np.double, np.single], np.single,replace=True)

