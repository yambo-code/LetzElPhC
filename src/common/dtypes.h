/*
This file contains the necessary data-structures
used in the Code.

Some comments :
1) "DO NOT" typedef a struct/enum. always use
    "struct"/"enum" keyword when defining
    for example
    correct way to define :
    ```
    struct lattice a ;
    struct wfc a ;
    enum ELPH_dft_code a;
    ```
*/
#pragma once
#include <mpi.h>
#include <stdbool.h>

#include "../elphC.h"
#include "error.h"

enum ELPH_dft_code
{
    DFT_CODE_QE
};

enum ELPH_screening
{
    ELPH_NO_SCREENING,
    ELPH_DFPT_SCREENING
};

struct kernel_info
{
    char name_str[32];
    bool bare_loc;
    bool non_loc;
    enum ELPH_screening screening;
};

struct symmetry
{
    ELPH_float Rmat[9];
    // rotation matrix
    ELPH_float tau[3];
    // fractional translation
    ND_int inv_idx;
    // index of inverse in the array.
    bool time_rev;
    // if this symmetry is time reversal
};

struct Lattice
{
    // type which contains all information about the crystal structure
    // -------------------------------------------------------
    int natom;
    // number of atoms in the unit cell
    int nsym;
    /* number of symmetries in nscf calculation (note this differs from
    symmetries used in phonon structure ) */
    int timerev;
    // System has time reversal symmetries
    int nspin;
    // number of spin components
    int nspinor;
    // number of spinor components
    int total_bands;
    // total number of bands used in nscf
    int start_band;
    // starting band used in elph calculation
    int end_band;
    // end band used in elph calculation
    int nbnds;
    // bands used in elph calculation
    char dimension;
    // '1', '2', '3' for 1D, 2D, 3D. This is used to apply coloumb cutoff
    ND_int nmag;
    /*  nmag = 1 for non magnetic case, nmag = 2 for lsda, and
    nmag = 4  for non-collinear magnetic */
    ND_int nmodes;
    // number of phonon modes
    ND_int nkpts_iBZ;
    // number of kpoints in iBZ
    ND_int nkpts_BZ;
    // number of kpoints in full BZ
    ND_int npw_max;
    // max number of total pw in spherical grid for all wfcs
    ND_int fft_dims[3];
    // fft dimensions
    ND_int nfftz_loc;
    // number of FFT vecs in this cpu
    ND_int nfftz_loc_shift;
    // global index of first fft vector
    ELPH_float alat_vec[9];
    // Lattice vectors a[i] = alat[:,i]
    ELPH_float blat_vec[9];
    // reciprocal lattice vectors b[i] = blat_vec[:,i] (includes 2*pi)
    ELPH_float volume;
    // Unit cell volume i.e det(alat_vec)
    ELPH_float* atomic_pos;
    // position of atoms in cartisian coordinates (natom,3)
    int* atom_type;
    // type of atoms array
    ELPH_float* kpt_iredBZ;
    // kpoints (cart coordinates) in irreducible BZ
    ELPH_float* kpt_fullBZ;
    // kpoints (cart coordinates) in full BZ
    ELPH_float* kpt_fullBZ_crys;
    // kpoints (crystal coordinates) in full BZ
    int* kmap;
    // map which relates the full BZ kpts to iBZ (nBZ,2), (kpt_iBZ,sym)
    struct symmetry* syms;
    // array of symmetries (operated on cartisian coordinates)
    bool is_soc_present;
    // system has soc if true else false
    // -------------------------------------------------------
};

struct Phonon
{
    // Struct to store some phonon related data
    // -------------------------------------------------------
    ND_int nq_iBZ;
    // number of points in iBZ
    ND_int nq_BZ;
    // number of points in full BZ
    ND_int nq_iBZ_loc;
    // q points in this q-pool
    ND_int nq_shift;
    /* shift to get global q idx.
    [nq_shift, nq_shift+nq_iBZ_loc] are in this pool */
    ND_int nph_sym;
    // number of phonon symmetries
    ELPH_float* qpts_iBZ;
    // qpoints in iBZ (crystal units)
    ELPH_float* qpts_BZ;
    // in full BZ (crystal units)
    struct symmetry* ph_syms;
    // Phonon symmetry matrices (cart)
    int* qmap;
    // (nq_BZ,2), (qpt_iBZ,nph_sym)
    ND_int* nqstar;
    // number of elements in each qstar of iBZ qpoints
    // -------------------------------------------------------
};

struct Vloc_table
{
    // Interpolation table for local part of the pseudo potential
    // in G space
    //  -------------------------------------------------------
    ELPH_float qmax_abs;
    /* absoule max bound of qpt coordinate (crystal coordinates)
    i.e, max(|q_{x,y,z}|) */
    ND_int npts_co;
    // number of points in the interpolation table
    ELPH_float* g_co;
    // |G| points (coarse grid)
    ELPH_float dg;
    // spacing between two g points in g_co
    ELPH_float* vlocg;
    // V(G) points on coarse grid (ntype,npts_co)
    ELPH_float* vploc_co;
    // V'(G) points on coarse grid (1st derivative) (ntype,npts_co)
    // -------------------------------------------------------
};

struct local_pseudo
{
    // contains data about local part of pseudo pots
    ELPH_float* Vloc_atomic;
    // local part of pseudo (ntype)
    ELPH_float* r_grid;
    // radial grid (ntype)
    ELPH_float* rab_grid;
    // derivative of radial grid = (ntype)
    ND_int ngrid;
    // number of points in the grid
    ELPH_float Zval;
    // number of valance electrons
};

struct Pseudo
{
    /* this structure contains all details about the pseudo potentials*/
    // -------------------------------------------------------
    // local data
    struct local_pseudo* loc_pseudo;
    // info about local parts (ntype)
    ELPH_float* PP_table;
    // PP_table from yambo. (nltimesj,natom_types,3) (l+1,2j,pp_spin), pp_spin =
    // 1 ?
    ELPH_float* Fsign;
    // sign of KB coeffients stored in "PP_KBS" (nlj_max, natom_types)
    ELPH_cmplx** fCoeff;
    // Refer Eq.9 of PHYSICAL REVIEW B 71, 115106 s2005d
    //  (natom_types,l*j) length array of (spin, spin, 2l+1, 2l+1)
    struct Vloc_table vloc_table[1];
    // interpolation table for local potential
    ND_int nltimesj;
    // nltimesj : max number of projectors (n*j)
    ND_int ngrid_max;
    // max value of grid of max(len(r_grid))
    ND_int ntype;
    // number of atom types
    int lmax;
    // Max value of l in the PP table
    // -------------------------------------------------------
};

/*struct for wfc */
struct WFC
{
    /* This struct is only used to store wfcs in iBZ */
    // -------------------------------------------------------
    ELPH_cmplx* wfc;
    // reciprocal space wavefunction (nspin,nbnd,nspinor,npw_loc)
    ELPH_float* gvec;
    // gvectors for which wfc is defined. (npw_loc,3)
    ELPH_float* Fk;
    // (nltimesj, ntype, npw_loc)
    // klein-bylander in reciprocal space (Fk in Eq. 6 of
    // https://docs.abinit.org/theory/pseudopotentials/)
    ND_int npw_total;
    // total number of gvecs for this wfc
    ND_int npw_loc;
    // gvecs in this cpu
    // -------------------------------------------------------
};

struct ELPH_MPI_Comms
{
    /* This struct contains all the information about communications
    used in the code*/
    MPI_Comm commW;
    // always set to MPI_COMM_"W"ORLD (do not free!)
    MPI_Comm commQ;
    // COMM for "Q"pool
    MPI_Comm commK;
    // COMM for "K"pool
    MPI_Comm commR;
    // equal "R"ank cpus from all kpools in COMM_WORLD form commR
    MPI_Comm commRq;
    // equal "R"ank cpus from all kpools in comm"Q" form commRk
    int nqpools;  // total number of q pools
    int nkpools;  // total number of k pools
    // ranks of comms
    int commW_rank;
    int commQ_rank;
    int commK_rank;
    int commR_rank;
    int commRq_rank;
    // size of comm
    int commW_size;
    int commQ_size;
    int commK_size;
    int commR_size;
    int commRq_size;
    /*
    for Example if there are total 8 cpus in MPI_COMM_WORLD
    divide as 2 qpools and 2 kpools

                            Global {0,1,2,3,4,5,6,7}
            _________________|__________________
            |                                   |
    {0,1,2,3} (commQ)                   {4,5,6,7} (commQ)
    _________|__________                _________|__________
    |                   |               |                   |
    {0,1}(commK)     {2,3}(commK)     {4,5} (commK)       {6,7} (commK)

    npw_cpus = size of commK
    commR contains  {0,2,4,6} and {1,3,5,7} i.e rank%npw_cpus is used to make
    the division commRq contains {0,2} , {1,3}, in 1st commQ and {4,5} and {5,7}
    in 2nd commQ
    */
};

/* This struct contains all the input file details */
struct usr_input
{
    // system varibles
    int nkpool;         // k point parallelization
    int nqpool;         // q point parallelization
    int start_bnd;      // starting band
    int end_bnd;        // last band
    char* save_dir;     // save dir
    char* ph_save_dir;  // ph_save directory
    char* kernel_str;   // level of screening to include
    bool kminusq;       // true if convention is "yambo" else false
};