/*
This file contains the necessary data-structures 
used in the Code.

Some comments :
1) "DO NOT" typedef a struct. always define the 
structure before prefixing "struct" for example
    write way to define :
    ```
    struct lattice a ;
    struct wfc a ;
    ```
2) structures should hold either pointers or a 
    single type. do not put array in structures
*/
#pragma once
#include "../elphC.h"

struct Lattice
{
    /* type which contains all information about the crystal structure 
    ---------------------------------  */
    int timerev  ;                          // System has time reversal symmetries
    int nspin ;                             // number of spin components
    int nspinor;                            // number of spinor components
    int total_bands ;                       // total number of bands used in nscf
    int start_band ;                        // starting band used in elph calculation
    int end_band ;                          // end band used in elph calculation
    int nbnds ;                             // bands used in elph calculation
    char dimension ;                        // '1', '2', '3' for 1D, 2D, 3D. This is used to apply coloumb cutoff
    ND_int npw_max;                         // max number of total pw in spherical grid for all wfcs
    ND_int fft_dims[3];                     // fft dimensions
    ND_int nffts_loc;                       // number of FFT vecs in this cpu
    ND_int nfft_shift_loc;                  // global index of first fft vector
    ND_array(Nd_floatS) * alat_vec ;        // Lattice vectors a[i] = alat[:,i]
    ND_array(Nd_floatS) * atomic_pos;       // position of atoms in cartisian coordinates
    int * atom_type ;                       // type of atoms array
    ND_array(Nd_floatS) * kpt_iredBZ ;      // kpoints (cart coordinates) in irreducible BZ
    ND_array(Nd_floatS) * kpt_fullBZ ;      // kpoints (cart coordinates) in full BZ
    ND_array(Nd_floatS) * kpt_fullBZ_crys ; // kpoints (crystal coordinates) in full BZ
    ND_array(i) * kmap ;                    // map which relates the full BZ kpts to iBZ (nBZ,2), (kpt_iBZ,sym)
    ND_array(Nd_floatS) * sym_mat ;         // symmetric matrices in cartesian system
    bool * time_rev_array ;                 // (nsym) bool array. True if the sym operation is Time rev
    ND_array(Nd_floatS) * frac_trans ;      /* fractional translation corresponding 
    ---------------------------------     to symmetric matrices (in cartesian coordinates) */

};

/* this structure contains all details about the pseudo potentials*/
struct Pseudo
{   
    // local data
    ND_array(Nd_floatS)* Vloc_atomic ; // local part of pseudo (ntype)

    ND_array(Nd_floatS)* r_grid ; //radial grid (ntype)

    ND_array(Nd_floatS)* rab_grid ; // derivative of radial grid = (ntype)

    // non local data
    ND_array(Nd_floatS) * PP_table ; // PP_table from yambo. (nlj_max,natom_types,3) (l+1,2j,pp_spin), pp_spin = 1 ?
    
    ND_array(Nd_floatS)* Fsign ; // sign of KB coeffients stored in "PP_KBS" (nlj_max,natom_types)
    
    ND_array(Nd_cmplxS) * fCoeff ; // Refer Eq.9 of PHYSICAL REVIEW B 71, 115106 s2005d 
    //  Coeff_len array of (spin, spin, 2l+1, 2l+1)

    ND_int ngrid_max ; // max value of grid of max(len(r_grid))

    ND_int ntype;      // number of atom types 

    ELPH_float * Zval; // Valance electrons for each atomic type

    int lmax;          // Max value of l in the PP table

};


/*struct for wfc */
struct WFC
{   /* This struct is only used to store wfcs in iBZ */
    ND_array(Nd_cmplxS) * wfc ;  
    // reciprocal space wavefunction (nspin,nbnd,nspinor,ng)
    
    ND_array(Nd_floatS) * gvec ; 
    // gvectors for which wfc is defined.

    ND_array(Nd_floatS) * Fk ;
    // klein-bylander in reciprocal space (Fk in Eq. 6 of https://docs.abinit.org/theory/pseudopotentials/)

    ND_int npw_total; // total number of gvecs for this wfc

    ND_int npw_loc;  // gvecs in this cpu

};



/*struct for fft_buffers */
struct wfcBox
{   
    /*
    buffer storage for wavefunctions in full cubic FFT grid and fft_plane for the buffers
    The rational behind this struct is that we create a optimized plan once 
    and reuse them multiple ttimes
    */
    ND_array(Nd_cmplxS) Buffer ; 
    // input buffer for FFT transform

    ND_array(Nd_cmplxS) Buffer_temp ; 
    // This is a temp buffer with size same as Buffer. used in local part calculation

    // this is a buffer to locally perform the FFT
    ND_array(Nd_cmplxS) FFTBuf ; // (nsets,Nx,Ny,Nz)

    // this is a buffer to store local wfc in Gsphere
    ELPH_cmplx * BufGsphere ; // (nsets,npw_total_max)

    ND_function(FFT_plan, Nd_cmplxS) fft_plan;
    // fft plan for the above arrays for forward transform

    ND_function(FFT_plan, Nd_cmplxS) inVfft_plan;
    //fft plan for the above arrays for forward transform

    ELPH_cmplx norm ; // normalization factor. This is just product of dims

    ELPH_float * Gvecs ; // temperary buffer to store Gvectors // (npw_total_max,3)

    ELPH_float * Gvecs_loc ; // local Gvecs buffers (npw_total_max/Comm + 1, 3)
    
    int * comm_buffer; // 4*Comm_size
    /*
    The use of the above two Gvec buffers and comm_buffer is we avoid repeative malloc calls 
    */
};
