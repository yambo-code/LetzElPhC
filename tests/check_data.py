## This file contains functions that checks the Data produced by Letzelphc

import numpy as np
from netCDF4 import Dataset
try:
    from scipy.spatial import cKDTree
except:
    from scipy.spatial import KDTree


rtol_elph = 1e-05  # relative tolerence used in allclose
atol_elph = 1e-8  # absoulte tolerence used in allclose

def set_tolerence(rtol=1e-05,atol=1e-8):
    rtol_elph = rtol
    atol_elph = atol


def quick_check_numeric_db(nc_db_test, nc_db_ref, var_name):
    # if return_data is true, then funtcion will test data (not the reference data)
    try:
        test_data = nc_db_test[var_name][...].data
        ref_data  = nc_db_ref[var_name][...].data
    except:
        return False
    return np.allclose(test_data,ref_data,rtol=rtol_elph, atol=atol_elph)
    

def make_kpositive(klist, tol=1e-6):
    ## brings all kpoints in [0,1)
    kpos = klist-np.floor(klist)
    return (kpos+tol)%1

#tree = spatial.KDTree(klist)
def find_kpt(kpt_search, tree, tol=1e-5):
    kpt_search = make_kpositive(kpt_search)
    dist, idx = tree.query(kpt_search, workers=1)
    if len(dist[dist>tol]) !=0:
        idx = -1
    return idx

def build_ktree(kpts):
    tree = make_kpositive(kpts)
    return cKDTree(tree,boxsize=[1,1,1])


def convert_yambo_to_std(eph_mat,kpts,qpts):
    ktree = build_ktree(kpts)
    # nq, nk, nmodes, nspin, initial_band, final_band_PH_ab
    nq, nk, nmodes = eph_mat.shape[:3]
    for iq in range(nq):
        idx_q = find_kpt(qpts[iq][None,:]+kpts, ktree)
        for imode in range(nmodes):
            elph_tmp_iq = eph_mat[iq,imode][idx_q,...].copy()
            eph_mat[iq,imode][...] = elph_tmp_iq[...]


def get_nc_strings(char_in):
    # Converts netcdf array of chars as python string.
    convlist=[iconv.decode('utf-8') for iconv in char_in]
    conv = ''
    for iconv in convlist: conv = conv+iconv
    return conv.strip()

def quick_check_char_db(nc_db_test, nc_db_ref, var_name):
    # if return_data is true, then funtcion will test data (not the reference data)
    try:
        test_data = nc_db_test[var_name][...].data
        ref_data  = nc_db_ref[var_name][...].data
    except:
        return False
    return (get_nc_strings(test_data) == get_nc_strings(ref_data))
    



def check_pol_vecs(nc_db_test, nc_db_ref):
    # We cannot directly compare the polarization vectors as they contain 
    # arbitary phases, instead, we construct the dynamical matrix and 
    # compare the dynamical matrices.
    try:
        test_data = nc_db_test['POLARIZATION_VECTORS'][...].data
        ref_data  = nc_db_ref['POLARIZATION_VECTORS'][...].data
        test_phfreq = nc_db_test['FREQ'][...].data
        ref_phfreq  = nc_db_ref['FREQ'][...].data
    except:
        return False
    
    test_data = test_data[...,0] + 1j*test_data[...,1]
    ref_data  = ref_data[...,0]  + 1j*ref_data[...,1]
        
    test_ph2 = np.where(test_phfreq < 0, -test_phfreq**2, test_phfreq**2)
    ref_ph2  = np.where(ref_phfreq < 0, -ref_phfreq**2, ref_phfreq**2)

    #
    if test_data.shape != ref_data.shape:
        return False
    
    nq, nmodes, atom, pol = ref_data.shape
    test_data = test_data.reshape(nq, nmodes, -1)
    ref_data = ref_data.reshape(nq, nmodes, -1)
    # nq, nmodes, atom, pol
    dyn_mat_ref = np.einsum('qv,qiv,qvj->qij',ref_ph2,np.linalg.inv(ref_data),ref_data,optimize=True)
    dyn_mat_test = np.einsum('qv,qiv,qvj->qij',test_ph2,np.linalg.inv(test_data),test_data,optimize=True)
    return np.allclose(dyn_mat_ref,dyn_mat_test,rtol=rtol_elph, atol=atol_elph)
    

def check_elph_me(nc_db_test, nc_db_ref):
    # We cannot directly compare the electron-phonon me as they contain 
    # arbitary phases from ph eigs vecs, so we need to multiply 
    # with overlap matrix between two phonon eigvecs
    try:
        test_qpts = nc_db_test['qpoints'][...].data
        ref_qpts  = nc_db_ref['qpoints'][...].data
        test_kpts = nc_db_test['kpoints'][...].data
        ref_kpts  = nc_db_ref['kpoints'][...].data
        test_polvec = nc_db_test['POLARIZATION_VECTORS'][...].data
        ref_polvec  = nc_db_ref['POLARIZATION_VECTORS'][...].data
        test_data = nc_db_test['elph_mat'][...].data
        ref_data  = nc_db_ref['elph_mat'][...].data
        ref_elph_convention  = nc_db_ref['convention'][...].data
        test_elph_convention = nc_db_test['convention'][...].data
    except:
        return False
    
    ref_elph_convention  = get_nc_strings(ref_elph_convention)
    test_elph_convention = get_nc_strings(test_elph_convention)

    test_data = test_data[...,0] + 1j*test_data[...,1]
    ref_data  = ref_data[...,0]  + 1j*ref_data[...,1]
    test_polvec = test_polvec[...,0] + 1j*test_polvec[...,1]
    ref_polvec  = ref_polvec[...,0]  + 1j*ref_polvec[...,1]
    if test_polvec.shape != ref_polvec.shape:
        return False
    if test_data.shape != ref_data.shape:
        return False
    nq, nmodes, atom, pol = ref_polvec.shape

    ## basic gauge invariance check
    test_data2 = np.sum(np.abs(test_data)**2,axis=(2,3,4,5))
    ref_data2  = np.sum(np.abs(ref_data)**2,axis=(2,3,4,5))

    test_pass0 = np.allclose(test_data2,ref_data2,rtol=rtol_elph)
    if not test_pass0:
        return test_pass0

    test_polvec = np.linalg.inv(test_polvec.reshape(nq, nmodes, -1))
    overlap_polvec = np.matmul(ref_polvec.reshape(nq, nmodes, -1),test_polvec)
    
    # incase if two conventions are not same, we need to convert to standard
    if (ref_elph_convention != test_elph_convention):
        if (ref_elph_convention != 'standard'):
            convert_yambo_to_std(ref_data, ref_kpts, ref_qpts)
        if (test_elph_convention != 'standard'):
            convert_yambo_to_std(test_data, test_kpts, test_qpts)
    # multiply with overlap
    test_data = np.einsum('qxv,qkvsij->qkxsij',overlap_polvec,test_data,optimize=True)
    ref_data_sum = np.sum(ref_data)
    
    test_pass1 = np.allclose(np.sum(test_data,axis=(3,4,5)),np.sum(ref_data,axis=(3,4,5)),rtol=rtol_elph, atol=1e-6)
    test_pass2 = test_pass1 and 100*abs(np.sum(test_data)-np.sum(ref_data))/abs(ref_data_sum) < 0.01

    return test_pass2 and np.allclose(test_data,ref_data,rtol=rtol_elph, atol=1e-6)


def check_elph_files(test_file, ref_file):
    test_pass = True
    try:
        elph_db_test = Dataset(test_file,'r')
        elph_db_ref  = Dataset(ref_file,'r')
    except:
        return False

    ## Check basic data
    test_pass = test_pass and quick_check_numeric_db(elph_db_test,elph_db_ref,'FREQ')        
    test_pass = test_pass and quick_check_numeric_db(elph_db_test,elph_db_ref,'kpoints')
    test_pass = test_pass and quick_check_numeric_db(elph_db_test,elph_db_ref,'qpoints')
    test_pass = test_pass and quick_check_numeric_db(elph_db_test,elph_db_ref,'qpoints_iBZ')
    test_pass = test_pass and quick_check_numeric_db(elph_db_test,elph_db_ref,'qmap')
    test_pass = test_pass and quick_check_numeric_db(elph_db_test,elph_db_ref,'kmap')
    test_pass = test_pass and quick_check_numeric_db(elph_db_test,elph_db_ref,'bands')
    test_pass = test_pass and quick_check_numeric_db(elph_db_test,elph_db_ref,'number_of_phonon_symmetries')
    test_pass = test_pass and quick_check_numeric_db(elph_db_test,elph_db_ref,'time_reversal_phonon')
    test_pass = test_pass and quick_check_numeric_db(elph_db_test,elph_db_ref,'symmetry_matrices')
    test_pass = test_pass and quick_check_numeric_db(elph_db_test,elph_db_ref,'fractional_translation')
    test_pass = test_pass and quick_check_char_db(elph_db_test,elph_db_ref, 'kernel')

    # check polarization vectors
    test_pass = test_pass and check_pol_vecs(elph_db_test,elph_db_ref)
    # check elph matrix elements
    test_pass = test_pass and check_elph_me(elph_db_test,elph_db_ref)

    elph_db_test.close()
    elph_db_ref.close()
    return test_pass


def check_dmat_files(test_file, ref_file):
    try:
        dmat_db_test = Dataset(test_file,'r')
        dmat_db_ref  = Dataset(ref_file,'r')
        test_data = dmat_db_test['Dmats'][...].data
        ref_data  = dmat_db_ref['Dmats'][...].data
        dmat_db_test.close()
        dmat_db_ref.close()
    except:
        return False
    
    test_data = test_data[...,0] + 1j*test_data[...,1]
    ref_data  = ref_data[...,0]  + 1j*ref_data[...,1]

    diff_max = np.abs(test_data-ref_data).max()
    return diff_max < 5e-4


# if (check_dmat_file("src/ndb.Dmats", "ndb.Dmats_ref_test")):
#     print("Test Passed")
# else :
#     print('Test Failed.')
