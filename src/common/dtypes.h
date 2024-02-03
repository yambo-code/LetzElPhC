/*
This file contains the necessary data-structures 
used in the Code.

Some comments :
1) "DO NOT" typedef a struct. always define the 
structure before prefixing "struct" for example
    correct way to define :
    ```
    struct lattice a ;
    struct wfc a ;
    ```
*/
#pragma once
#include "../elphC.h"

struct Lattice
{
    // type which contains all information about the crystal structure 
    // ------------------------------------------------------- 
    int timerev  ;                          
    // System has time reversal symmetries
    int nspin ;                             
    // number of spin components
    int nspinor;                            
    // number of spinor components
    int total_bands ;                       
    // total number of bands used in nscf
    int start_band ;                        
    // starting band used in elph calculation
    int end_band ;                          
    // end band used in elph calculation
    int nbnds ;                             
    // bands used in elph calculation
    char dimension ;                        
    // '1', '2', '3' for 1D, 2D, 3D. This is used to apply coloumb cutoff
    ND_int npw_max;                         
    // max number of total pw in spherical grid for all wfcs
    ND_int fft_dims[3];                     
    // fft dimensions
    ND_int nfftz_loc;                       
    // number of FFT vecs in this cpu
    ND_int nfftz_loc_shift;                 
    // global index of first fft vector
    ND_array(Nd_floatS) * alat_vec ;        
    // Lattice vectors a[i] = alat[:,i]
    ND_array(Nd_floatS) * atomic_pos;       
    // position of atoms in cartisian coordinates
    int * atom_type ;                       
    // type of atoms array
    ND_array(Nd_floatS) * kpt_iredBZ ;      
    // kpoints (cart coordinates) in irreducible BZ
    ND_array(Nd_floatS) * kpt_fullBZ ;      
    // kpoints (cart coordinates) in full BZ
    ND_array(Nd_floatS) * kpt_fullBZ_crys ; 
    // kpoints (crystal coordinates) in full BZ
    ND_array(i) * kmap ;                    
    // map which relates the full BZ kpts to iBZ (nBZ,2), (kpt_iBZ,sym)
    ND_array(Nd_floatS) * sym_mat ;         
    // symmetric matrices in cartesian system
    bool * time_rev_array ;                 
    // (nsym) bool array. True if the sym operation is Time rev
    ND_array(Nd_floatS) * frac_trans ;      
    /* fractional translation corresponding 
    to symmetric matrices (in cartesian coordinates) */
    bool is_soc_present;                    
    // system has soc if true else false
    ND_int nmag;                            
    /*  nmag = 1 for non magnetic case, nmag = 2 for lsda, and
    nmag = 4  for non-collinear magnetic */
    // ------------------------------------------------------- 
};

struct Phonon
{   
    //      Struct to store some phonon related data
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
    ELPH_float * qpts_iBZ;      
    // qpoints in iBZ (crystal units)
    ELPH_float * qpts_BZ;       
    // in full BZ (crystal units)
    ELPH_float * ph_sym_mats;   
    // Phonon symmetry matrices
    ELPH_float * ph_sym_tau;    
    // Corresponding frac. trac. vecs
    bool       * time_rev_array ; 
    // (nsym) bool array. True if the sym operation is Time rev
    int * qmap ;
    // (nq_BZ,2), (qpt_iBZ,nph_sym)
    ND_int * nqstar; 
    // number of elements in each qstar of iBZ qpoints
    // ------------------------------------------------------- 
};


struct Vloc_table
{
    //Interpolation table for local part of the pseudo potential 
    //in G space
    // ------------------------------------------------------- 
    ELPH_float qmax_abs;
    /* absoule max bound of qpt coordinate (crystal coordinates) 
    i.e, max(|q_{x,y,z}|) */
    ND_int npts_co; 
    // number of points in the interpolation table
    ELPH_float * g_co;  
    // |G| points (coarse grid)
    ELPH_float dg; 
    // spacing between two g points in g_co
    ELPH_float * vlocg; 
    // V(G) points on coarse grid (ntype,npts_co)
    ELPH_float * vploc_co; 
    // V'(G) points on coarse grid (1st derivative) (ntype,npts_co)
    // ------------------------------------------------------- 
};


struct Pseudo
{   
    /* this structure contains all details about the pseudo potentials*/
    // ------------------------------------------------------- 
    // local data
    ND_array(Nd_floatS)* Vloc_atomic ; 
    // local part of pseudo (ntype)
    ND_array(Nd_floatS)* r_grid ; 
    //radial grid (ntype)
    ND_array(Nd_floatS)* rab_grid ; 
    // derivative of radial grid = (ntype)
    // non local data
    ND_array(Nd_floatS) * PP_table ; 
    // PP_table from yambo. (nlj_max,natom_types,3) (l+1,2j,pp_spin), pp_spin = 1 ?
    ND_array(Nd_floatS)* Fsign ; 
    // sign of KB coeffients stored in "PP_KBS" (nlj_max,natom_types)
    ND_array(Nd_cmplxS) * fCoeff ; 
    // Refer Eq.9 of PHYSICAL REVIEW B 71, 115106 s2005d 
    //  Coeff_len array of (spin, spin, 2l+1, 2l+1)
    struct Vloc_table vloc_table[1]; 
    // interpolation table for local potential
    ND_int ngrid_max ; 
    // max value of grid of max(len(r_grid))
    ND_int ntype;      
    // number of atom types 
    ELPH_float * Zval; 
    // Valance electrons for each atomic type
    int lmax;          
    // Max value of l in the PP table
    // ------------------------------------------------------- 
};


/*struct for wfc */
struct WFC
{   
    /* This struct is only used to store wfcs in iBZ */
    // ------------------------------------------------------- 
    ND_array(Nd_cmplxS) * wfc ;  
    // reciprocal space wavefunction (nspin,nbnd,nspinor,ng)
    ND_array(Nd_floatS) * gvec ; 
    // gvectors for which wfc is defined.
    ND_array(Nd_floatS) * Fk ;
    // klein-bylander in reciprocal space (Fk in Eq. 6 of https://docs.abinit.org/theory/pseudopotentials/)
    ND_int npw_total; // total number of gvecs for this wfc
    ND_int npw_loc;  // gvecs in this cpu
    // ------------------------------------------------------- 
};


