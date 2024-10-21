from netCDF4 import Dataset
import numpy as np
from matplotlib import pyplot as plt
from numpy import einsum
from multiprocessing import Pool
import itertools
import glob
import sys 
np.set_printoptions(suppress=True)


path_dir = sys.argv[1]
code = sys.argv[2]
qe = False
if code.strip() == 'qe':
    qe = True

path_dir = path_dir.strip()
prefix = glob.glob(path_dir+"/*.save")[0]
ns_db = prefix.strip()+'/SAVE/ns.db1'
dipole_file = prefix.strip()+'/SAVE/ndb.dipoles'

if qe:
    elph_file = path_dir+'/elph_dir/elph_yambo_s.dbph_000001.h5'
else :
    elph_file = path_dir+'/ndb.elph'


############ kpoint expansion stuff #######

def k_is_present(kpoints, kpt, tol=10**-4):
    if len(kpoints) == 0:
        return False
    k_temp_points = np.array(kpoints) - np.array(kpt)[None, :]
    k_temp_points = k_temp_points-np.rint(k_temp_points)
    k_temp_points = np.linalg.norm(k_temp_points, axis=-1)
    return np.isclose(k_temp_points, 0, atol=tol).any()
# expand kpoints to full brillioun zone

def expand_kpoints(kibz, symmetries, lat_vec):
    # kibz is cart corrdinates
    ## This also brings kpoints to wigner-sieze cell.
    #S^T@k.T (nbz,3)  k@S
    b_vec = np.linalg.inv(lat_vec)
    k_full = np.einsum('kj,nji->nki', kibz, symmetries, optimize=True)
    kpoints = []
    kmap = []
    temp_range = np.arange(-2, 3)
    combinations = np.array(list(itertools.product(temp_range, temp_range, temp_range)))@b_vec  # cart coordinates
    for i in range(len(kibz)):
        for j in range(len(symmetries)):
            temp = np.linalg.norm(k_full[j, i][None, :]-combinations, axis=-1).argmin()
            k_full[j, i] = k_full[j, i]-combinations[temp]
            k_cyrs = k_full[j, i][None, :]@lat_vec
            k_cyrs = k_cyrs.reshape(-1)
            if len(kpoints) > 0:
                if not k_is_present(kpoints, k_cyrs):
                    kpoints.append(k_cyrs)
                    kmap.append([i, j])
            else:
                kpoints.append(k_cyrs)
                kmap.append([i, j])
    return np.array(kpoints)@b_vec, np.array(kmap)




##########  reading stuff


db = Dataset(ns_db, mode='r')
dims_db = db['DIMENSIONS'][...].data 
lat_param = db['LATTICE_PARAMETER'][...].data
lat_vec = db['LATTICE_VECTORS'][...].data
yambo_kpts = db['K-POINTS'][...].data.T
symmetries = db['SYMMETRY'][...]
yambo_kpts = yambo_kpts/lat_param[None, :] ## yambo_kpts in cart
kpts_BZ, kmap = expand_kpoints(yambo_kpts, symmetries, lat_vec)
time_rev = int(np.rint(dims_db[9]))
nsym = symmetries.shape[0]
mid_sym = nsym//(time_rev+1)
time_rev_idx = (kmap[:,1] >= mid_sym)


dipoles = Dataset(dipole_file)['DIP_v'][...][0]
# (nspin, nk, nv, nc ,npol, re/im)
Nk     = dipoles.shape[0]
Nv_bse = dipoles.shape[1]
Nc_bse = dipoles.shape[2]

dipoles = dipoles[...,0] - 1j*dipoles[...,1] # // nk, nv, nc ,npol
dipoles = dipoles.transpose(3,0,2,1) #// npol, nk, nc, nv

dip_v = dipoles

pars = Dataset(dipole_file)['PARS'][...].data

nvmax =  int(pars[2])
dip_v_full_bz = dip_v[:,kmap[:,0],...]
dip_v_full_bz[:,time_rev_idx] = np.conj(dip_v_full_bz[:,time_rev_idx])
sym_mats = symmetries[kmap[:,1],...]
dip_v = np.einsum('kji,jkcv->ikcv',sym_mats,dip_v_full_bz,optimize=True)
Nk = dip_v.shape[1]

#print(kpts_BZ@lat_vec)
print('Total kpoints in BZ : ', Nk)

elph = Dataset(elph_file)
if qe:
    elph = elph['ELPH_GKKP_Q1'][...].data
    elph = elph[...,0] + 1j*elph[...,1]
    elph = elph.transpose(3,0,2,1)
else :
    elph = elph['elph_mat'][...].data*(0.5)**1.5
    elph = elph[...,0] + 1j*elph[...,1] #// "nk","nmodes","nspin","mk","nkq" 
    elph = elph[0,:,:,0,:,:]
    elph = elph.transpose(1,0,3,2) # // (branch) (k-point) (final band) (initial band)


g = elph[:,:,:,:]
nmodes = g.shape[0]
print(Nv_bse)

d_mu = dip_v

sum = np.einsum('xkvc,ikcv->xi',g[:,:,:Nv_bse,Nv_bse:],d_mu)
# gcc = g[:,:,Nv_bse:,Nv_bse:]
# gvv = g[:,:,:Nv_bse,:Nv_bse]

# g_temp = einsum('xkci,nkiv->xnkcv',gcc,d_mu) #// c'c, cv->c'v
# g_temp = np.conj(g_temp) + einsum('xkvi,nkci->xnkcv',gvv,np.conj(d_mu)) # vv', cv->v'c


# sum = einsum('mkcv,xnkcv->xmn',d_mu,g_temp)
print(np.sum(np.abs(sum)**2))  ### sum is gauge invariant
print(np.sum(np.abs(g)**2,axis=(1,2,3)))

