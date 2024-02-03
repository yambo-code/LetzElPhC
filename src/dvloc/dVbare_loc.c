#include "dvloc.h"


void dVlocq(const ELPH_float * qpt, struct Lattice * lattice, struct Pseudo * pseudo, \
            ND_array(Nd_cmplxS) * eigVec, ND_array(Nd_cmplxS) * Vlocr, MPI_Comm commK)
{   
    /* 
    Computes the change in local bare potential in real space
    
    Input 
    qpt       : q-point in crystal coordinates (FIX ME. qpt must be in [-1,1]
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
    
    int mpi_error;
    // get the total cpus in comm and rank of each cpu
    int my_rank, Comm_size;
    mpi_error = MPI_Comm_size(commK, &Comm_size);
    mpi_error = MPI_Comm_rank(commK, &my_rank);


    const ND_array(Nd_floatS)* latvec       = lattice->alat_vec;
    const ND_int ntype                      = pseudo->ntype;
    const ND_array(Nd_floatS)* atom_pos     = lattice->atomic_pos;
    const int * atom_type                   = lattice->atom_type; 
    const char cutoff                       = lattice->dimension;
    const ND_int ngrid_max                  = pseudo->ngrid_max;
    const ELPH_float * Zval                 = pseudo->Zval;
    ND_int natom                            = atom_pos->dims[0];

    ND_function(set_all, Nd_cmplxS)(Vlocr , 0.0);

    ND_int G_vecs_xy, Gxy_shift;
    
    G_vecs_xy = get_mpi_local_size_idx(lattice->fft_dims[0]*lattice->fft_dims[1], \
                                        &Gxy_shift, commK);

    ND_int size_G_vecs      = G_vecs_xy*lattice->fft_dims[2]; 
    ND_int size_VG          = 3*natom*size_G_vecs;

    ELPH_cmplx *  VlocG      = malloc(size_VG*sizeof(ELPH_cmplx)); // 3*natom* ix_s*jy_s*kz
    
    int * gvecs = malloc(3*size_G_vecs*sizeof(int));
    
    ELPH_float volume = fabs(det3x3(latvec->data));

    ELPH_float blat[9]; // 2*pi is included 
    reciprocal_vecs(latvec->data,blat);
    
    /* This is spread of nuclear charge */
    ELPH_float eta = 1.0 ;

    ND_int FFTx, FFTy, FFTz;
    FFTx = lattice->fft_dims[0];
    FFTy = lattice->fft_dims[1];
    FFTz = lattice->fft_dims[2];

    ND_int npts_co          = pseudo->vloc_table->npts_co; 
    ELPH_float * g_co       = pseudo->vloc_table->g_co;
    ELPH_float * vlocg      = pseudo->vloc_table->vlocg;
    ELPH_float * vploc_co   = pseudo->vloc_table->vploc_co;

    ELPH_float dg = pseudo->vloc_table->dg; // spacing between two g points in g_co

    ELPH_float * VlocGtype = malloc(sizeof(ELPH_float)*(ntype));
    
    
    for (ND_int ig = 0 ; ig < G_vecs_xy; ++ig)
    {   
        ND_int fft_glob_idx = Gxy_shift + ig;
        ND_int ix = fft_glob_idx/FFTy ;
        ND_int jy = fft_glob_idx%FFTy ;
        
        for (ND_int kz = 0; kz<FFTz ; ++kz )
        {   
            int * g_temp_set =  gvecs + FFTz*3*ig + kz*3;
            g_temp_set[0] = ix; g_temp_set[1] = jy;  g_temp_set[2] = kz;

            ELPH_float qGtemp[3] = {get_miller_idx(ix,FFTx)+qpt[0], \
            get_miller_idx(jy,FFTy)+qpt[1], get_miller_idx(kz,FFTz)+qpt[2]}; // | q + G|

            ELPH_float qGtempCart[3]; // in cartisian coordinate

            MatVec3f(blat,qGtemp,false,qGtempCart); // 2*pi is included here //
                
            ELPH_float qGnorm = sqrt(qGtempCart[0]*qGtempCart[0] + \
                                qGtempCart[1]*qGtempCart[1] + qGtempCart[2]*qGtempCart[2]);
            
            if (qGnorm > g_co[npts_co-1]) error_msg("Vlocg Interpolation failed due to out of bounds");

            ND_int gidx = floor(qGnorm/dg);

            ELPH_cmplx * tmp_ptr =  VlocG + 3*natom*(ig*FFTz + kz); 

            ELPH_float cutoff_fac = 1;
            /* using analytic cutoff which works only when z periodicity is broken */
            if (cutoff == '2')
            {   
                ELPH_float qGp = latvec->data[8] * sqrt(qGtempCart[0]*qGtempCart[0] + qGtempCart[1]*qGtempCart[1]) ; 
                cutoff_fac -= exp(-qGp*0.5)*cos(qGtempCart[2]*latvec->data[8]*0.5) ;
            }

            for (ND_int itype = 0; itype <ntype; ++itype )
            {
                VlocGtype[itype] = spline_interpolate(qGnorm, gidx, g_co, vlocg+itype*npts_co , vploc_co+itype*npts_co);
                // add long range part back
                if (qGnorm >= ELPH_EPS) 
                {
                    VlocGtype[itype] -= \
                    (4*ELPH_PI*ELPH_e2)*Zval[itype]*exp(-qGnorm*qGnorm*0.25/eta)*cutoff_fac/(qGnorm*qGnorm*volume) ;
                }
                //printf("%f \n",VlocGtype[itype]);
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
    }

    free(VlocGtype);

    ELPH_cmplx *  VlocG_mode = calloc(size_VG, sizeof(ELPH_cmplx)); // 3*natom* ix_s*jy_s*kz

    //VlocG -> (nffs,atom,3), 
    // get it in mode basis @ (nu,atom,3) @(nffs,atom,3)^T
    ND_int ldA_eigVec = eigVec->strides[0]; //
    ND_int nmodes_ph  = eigVec->dims[0];


    ND_function(matmulX, Nd_cmplxS) ('N', 'T', eigVec->data, VlocG, VlocG_mode, \
        1.0, 0.0, ldA_eigVec, ldA_eigVec, size_G_vecs, nmodes_ph, size_G_vecs, ldA_eigVec);
    free(VlocG);

    struct ELPH_fft_plan fft_plan;

    wfc_plan(&fft_plan, size_G_vecs, lattice->nfftz_loc, G_vecs_xy, gvecs, \
            lattice->fft_dims, FFTW_MEASURE, commK);

    invfft3D(&fft_plan, nmodes_ph, VlocG_mode, Vlocr->data, false);
    
    wfc_destroy_plan(&fft_plan);
    // ND_function(write, Nd_cmplxS) ("Vloc.nc","local_pot",Vlocr,(char * []){"mode","ix","jy","kz"},NULL ) ;
    free(gvecs);
    free(VlocG_mode);

    
}

