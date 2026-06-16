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
};

