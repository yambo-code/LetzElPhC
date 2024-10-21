# coding: utf-8
from numpy import *
import numpy as np


qpts = np.loadtxt('qpts.in')
kpts = np.loadtxt('kpts') 

from netCDF4 import Dataset

def get_nc(file,var,cmplx=True):
    aaaaaa = Dataset(file)[var][...].data
    if cmplx:
        return aaaaaa[...,0] + 1j*aaaaaa[...,1]
    else :
        return aaaaaa

kminusq_idx = []

for iq in range(len(qpts)):
    kmq = kpts-qpts[iq][None,:]
    kkkkk = kmq[:,None,:]- kpts[None,:,:]
    kkkkk = kkkkk-np.rint(kkkkk)
    kkkkk = np.linalg.norm(kkkkk,axis=-1)
    kminusq_idx.append(kkkkk.argmin(axis=-1))

kminusq_idx = np.array(kminusq_idx)

elph_yambo = get_nc("ndb.elph_yambo","elph_mat")
elph_stand = get_nc("ndb.elph_standard","elph_mat")
print(np.abs(elph_stand-elph_yambo).max())
elph_yambo_new = np.zeros(elph_yambo.shape,dtype=complex)
for iq in range(len(qpts)):
    elph_yambo_new[iq][kminusq_idx[iq]] = elph_yambo[iq]

#print(elph_yambo_new.shape)
#print(elph_stand[:, 14, 5, 0, 4, 5])
#print(elph_yambo_new[:, 14, 5, 0, 4, 5])
print(np.abs(elph_stand-elph_yambo_new).max())
