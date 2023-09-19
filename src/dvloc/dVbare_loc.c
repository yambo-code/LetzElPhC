#include "dvloc.h"

static ND_int get_miller_idx(ND_int idx_in, ND_int FFT_dimension);


void dVlocq(const ELPH_float * qpt, struct Lattice * lattice, struct Pseudo * pseudo, \
            ND_array(Nd_cmplxS) * eigVec, ND_array(Nd_cmplxS) * Vlocr, MPI_Comm commK)
{   
    /* 
    Computes the change in local bare potential on reciprocal fft grid. the output must 
    be fourier transformed
    
    Input 
    qpt       : q-point in crystal coordinates
    atom_pos  : atomic positions in cart coordinates
    nDim      : '3','2','1' for 3D, 2D, 1D cutoffs
    Vloc_atomic :   Local potential of eqch atom in radial grid.
                    Note this is read from pseudo potentials and is 
                    array of ND_arrays (ntype) arrays of (radial grid)
    rab_grid  : drab (read from pseudos). used for integration. (ntype)
    Zval      : (ntype) array of valance charge
    eigVec    : eigen vectors , (nu,atom,3)
    atom_type : atomic type array (natom) for ex: Zval[atom_type[i]] gives
                atomic number of atom i
    Out-put:
            Vlocr
            d Vloc/dtau (in cart) (nmode,nffts_this_cpu), 
    */

    const ND_array(Nd_floatS)* latvec       = lattice->alat_vec;
    const ND_int ntype                      = pseudo->ntype;
    const ND_array(Nd_floatS)* atom_pos     = lattice->atomic_pos;
    const int * atom_type                   = lattice->atom_type; 
    const char cutoff                       = lattice->dimension;
    const ND_array(Nd_floatS)* Vloc_atomic  = pseudo->Vloc_atomic;
    const ND_array(Nd_floatS)* r_grid       = pseudo->r_grid;
    const ND_array(Nd_floatS)* rab_grid     = pseudo->rab_grid;
    const ND_int ngrid_max                  = pseudo->ngrid_max;
    const ELPH_float * Zval                 = pseudo->Zval;

    ND_int size_Vr      = ND_function(size, Nd_cmplxS) (Vlocr);
    ELPH_cmplx *  VlocG = malloc(sizeof(ELPH_cmplx)*size_Vr);

    ND_function(set_all, Nd_cmplxS)(Vlocr , 0.0);
    for (ND_int ig = 0; ig < size_Vr; ++ig) VlocG[ig] = 0.0 ;
    
    ND_int natom = atom_pos->dims[0];
    
    ELPH_float volume = fabs(det3x3(latvec->data));

    ELPH_float blat[9]; // 2*pi is included 
    reciprocal_vecs(latvec->data,blat);
    
    /* This is spread of nuclear charge */
    ELPH_float eta = 1.0 ;

    ND_int FFTx, FFTy, FFTz;
    FFTx = lattice->fft_dims[0];
    FFTy = lattice->fft_dims[1];
    FFTz = lattice->fft_dims[2];

    ND_int fft_strides[3] = {FFTy*FFTz,FFTz,1};

    ELPH_OMP_PAR
    {
    ELPH_float * work_array;
    ELPH_OMP_PAR_CRITICAL
    work_array = malloc(sizeof(ELPH_float)*(ntype+ngrid_max));
    
    ELPH_float * VlocGtype = work_array+ngrid_max;
    
    for (ND_int ifft = 0 ; ifft < lattice->nffts_loc; ++ifft)
    {   
        ND_int fft_glob_idx = lattice->nfft_shift_loc + ifft;
        ND_int Nx = fft_glob_idx/fft_strides[0] ;
        ND_int temp_rem = fft_glob_idx%fft_strides[0] ;
        ND_int Ny = temp_rem/fft_strides[1] ;
        ND_int Nz = temp_rem%fft_strides[1] ;
        

        ELPH_float qGtemp[3] = {get_miller_idx(Nx,FFTx)+qpt[0], \
        get_miller_idx(Ny,FFTy)+qpt[1], get_miller_idx(Nz,FFTz)+qpt[2]}; // | q + G|

        ELPH_float qGtempCart[3]; // in cartisian coordinate

        MatVec3f(blat,qGtemp,false,qGtempCart); // 2*pi is included here //
                
        ELPH_float qGnorm = sqrt(qGtempCart[0]*qGtempCart[0] + \
                            qGtempCart[1]*qGtempCart[1] + qGtempCart[2]*qGtempCart[2]);
                            
        ELPH_cmplx * tmp_ptr =  VlocG + 3*natom*ifft; 

        ELPH_float cutoff_fac = 1;
        /* using analytic cutoff which works only when z periodicity is broken */
        if (cutoff == '2')
        {   
            ELPH_float qGp = latvec->data[8] * sqrt(qGtempCart[0]*qGtempCart[0] + qGtempCart[1]*qGtempCart[1]) ; 
            cutoff_fac -= exp(-qGp*0.5)*cos(qGtempCart[2]*latvec->data[8]*0.5) ;
        }

        for (ND_int itype = 0; itype <ntype; ++itype )
        {
            VlocGtype[itype] = Vloc_Gspace(work_array, cutoff, qGnorm, Vloc_atomic+itype, \
                            r_grid+itype, rab_grid+itype, Zval[itype],eta,cutoff_fac)/volume ;
        }
        ELPH_OMP_PAR_SIMD
        for (ND_int ia = 0 ; ia<natom; ++ia)
        {   
            ND_int itype = atom_type[ia];
                    
            ELPH_float * pos_temp = atom_pos->data + 3*ia ;
            ELPH_float qdottau = qGtempCart[0]*pos_temp[0]+qGtempCart[1]*pos_temp[1] + qGtempCart[2]*pos_temp[2];
            ELPH_cmplx factor = -I*cexp(-I*qdottau);
            factor *= VlocGtype[itype];
            tmp_ptr[ia*3]   = factor*qGtempCart[0];
            tmp_ptr[ia*3+1] = factor*qGtempCart[1];
            tmp_ptr[ia*3+2] = factor*qGtempCart[2];
        }
    }
    ELPH_OMP_PAR_CRITICAL
    free(work_array);
    }
    //VlocG -> (nffs,atom,3), 
    // get it in mode basis @ (nu,atom,3) @(nffs,atom,3)^T
    ND_int ldA_eigVec = eigVec->strides[0]; //
    ND_int nmodes_ph  = eigVec->dims[0];
    ND_int nfft_dim   = lattice->nffts_loc ;

    ND_function(matmulX, Nd_cmplxS) ('N', 'T', eigVec->data, VlocG, Vlocr->data, \
        1.0, 0.0, ldA_eigVec, ldA_eigVec, nfft_dim, nmodes_ph, nfft_dim, ldA_eigVec);
    
    // free some buffers
    free(VlocG);
    
    // perform the inv FFT
    VlocinVFFT(Vlocr, lattice, commK);
    
}




/* Static helper functions */
static ND_int get_miller_idx(ND_int idx_in, ND_int FFT_dimension)
{
    // returns FFT Indices to [-N/2, N/2) if FFT_dimension is even
    // [-(n-1)/2,(n-1)/2 )] of FFT_dimension is odd
    ND_int mid_pnt = (FFT_dimension-1)/2 + 1 ;
    if ( idx_in < mid_pnt ) return idx_in;
    else return (idx_in-mid_pnt)-FFT_dimension/2 ;
}