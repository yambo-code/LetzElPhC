/*
This file contains Yambo specific data-structures used in the Code.
*/
#pragma once
#include <mpi.h>

struct Y6_parallel_work
{
 int *Q;
 int NQ;
 int *K;
 int NK;
};

// Pieces of Lattice and Phonon to be passed to Yambo//
struct Y6_info
{
    int natom;
    int nsym;
    int timerev;
    int nspin;
    int nspinor;
    int total_bands;
    int start_band;
    int end_band;
    int nbnds;
    int lattice_dim;
    int nmag;
    int nmodes;
    int nq_ibz;
    int nq_bz;
    int nkpts_ibz;
    int nkpts_bz;
    int nph_sym;
    int kminusq;
    // Pointers to k/q lists (C-side allocated, valid for query lifetime)
    float* kpt_fullBZ_crys;  // (nkpts_bz, 3) full BZ k-points in crystal coords
    int* kmap;               // (nkpts_bz, 2) maps full BZ k to iBZ+sym
    float* qpts_iBZ;         // (nq_ibz, 3) iBZ q-points
    float* qpts_BZ;          // (nq_bz, 3) full BZ q-points
    int* qmap;               // (nq_bz, 2) maps full BZ q to iBZ+sym
    int* nqstar;             // (nq_ibz,) star multiplicity
    int* kplusq_idxs;        // (nq_ibz, nkpts_bz) k+q indices for each q,k pair
};

