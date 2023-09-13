#include "dvloc.h"


void compute_elphLocal_q(struct WFC * wfcs, struct Lattice * lattice, struct Pseudo * pseudo, \
            ELPH_float * qpt, ND_array(Nd_cmplxS) * eigVec, ND_array(Nd_cmplxS) * dVscf_full, 
            ELPH_cmplx * elph_kq)
{
    /*
    dVscf -> (nmodes,nmag,Nx,Ny,Nz)
    ((k, nmodes, nspin, nbands, nbands))
    */

    ND_int * FFT_dims = dVscf_full->dims+2;
    ND_int iqpt ;
    
    bool qis_present = Function(isVECpresent, Nd_cmplxS) (lattice->kpt_fullBZ_crys,qpt,&iqpt);
    
    if (!qis_present) error_msg("q point is not commensurate with the kpoints");

    int q_ibz_idx  = *ND_function(ele,i) (lattice->kmap,nd_idx{iqpt,0});
    int q_ibz_sym  = *ND_function(ele,i) (lattice->kmap,nd_idx{iqpt,1});

    const ELPH_float * iqsym = ND_function(ele,Nd_floatS)(lattice->sym_mat, nd_idx{q_ibz_sym,0,0}) ;

    ND_int nmodes = dVscf_full->dims[0];

    ND_int ngvecs = (wfcs+q_ibz_idx)->gvec->dims[0];

    ND_array(Nd_cmplxS) VlocG, VlocG_iv, dVscf ;

    /** FIX ME need to removed when reading fron netcdf*/
    ND_function(init_strip_dims, Nd_cmplxS) (dVscf_full, 1, &dVscf);

    ND_function(init,Nd_cmplxS)(&VlocG, 2, nd_idx{nmodes,ngvecs});

    ND_function(init,Nd_cmplxS)(&VlocG_iv, 4, nd_idx{1,1,1,ngvecs});

    ND_function(malloc,Nd_cmplxS)(&VlocG); // internally, it will be zeroed in dVGlocq
    
    dVGlocq(wfcs, lattice, pseudo, iqpt, eigVec,&VlocG);

    const ELPH_float zero_temp[3] = {0.0,0.0,0.0} ;
    ND_int mag_iter = 1;
    if (lattice->nspin == 2) mag_iter = 2;

    int * kmap  = lattice->kmap->data ; 
    int * KplusQidxs = malloc((lattice->kmap->dims[0])*sizeof(int));
    get_KplusQ_idxs(lattice->kpt_fullBZ_crys, KplusQidxs , qpt, lattice->alat_vec, true);

    ND_int nbnd = wfcs->wfc->dims[1] ; 
    ND_int elph_kstride_mode = lattice->nspin *nbnd * nbnd; // 1st stride value of elph_kq
    ND_int elph_kstride_k    = eigVec->dims[0] *elph_kstride_mode;
    struct wfcBox wfcr_buffer;
    //printf("%lld %lld %lld \n",FFT_dims[0], FFT_dims[1], FFT_dims[2]);
    alloc_wfcBox(&wfcr_buffer,  6, nd_idx{1, 1, 1, FFT_dims[0], FFT_dims[1], FFT_dims[2]}, 3, nd_idx{3,4,5}, true, FFTW_MEASURE);
            
    for (ND_int iv = 0 ; iv < nmodes; ++iv)
    {   
        dVscf.data = dVscf_full->data + iv*dVscf_full->strides[0];

        VlocG_iv.data = VlocG.data + iv*ngvecs ;

        ND_int nmag = dVscf.dims[0];

        wfcinVFFT(&VlocG_iv, NULL, iqsym, zero_temp, NULL, false, \
        false, false, lattice->alat_vec->data, (wfcs+q_ibz_idx)->gvec, &wfcr_buffer);

        /* Note that we assume we do not have any external magnetic field, so electric field from nuclei*/
        for (ND_int im = 0; im<mag_iter; ++im)
        {
            ELPH_cmplx * restrict temp_ptr = dVscf.data + (im*dVscf.strides[0]) ; // mxyzai + xyzai->
            ELPH_OMP_PAR_FOR_SIMD
            for (ND_int i = 0 ; i< dVscf.strides[0]; ++i) temp_ptr[i] += wfcr_buffer.outBuf.data[i] ;
        }
        if(nmag == 4)
        {
            /* dvscf_{2x2} = vloc*I + Bx*sigma_x + By*sigma_y + Bz*sigma_z
            | V + Bz         Bx - I*By |
            | Bx+ I*By       V - Bz    |
            where Bx,By,Bz = d{Exc}/dm and Vxc = d{Exc}/dn
            */

            /*
            This is an inplace operation to avoid creating a large buffer
            */

            ELPH_cmplx * restrict Vxc_loc = dVscf.data + (0*dVscf.strides[0]) ;
            ELPH_cmplx * restrict  Bx_loc = dVscf.data + (1*dVscf.strides[0]) ;
            ELPH_cmplx * restrict  By_loc = dVscf.data + (2*dVscf.strides[0]) ;
            ELPH_cmplx * restrict  Bz_loc = dVscf.data + (3*dVscf.strides[0]) ;
            ELPH_OMP_PAR_FOR_SIMD
            for (ND_int i = 0 ; i< dVscf.strides[0]; ++i)
            {
                ELPH_cmplx Vxc_t, Bx_t, By_t, Bz_t;
                Vxc_t = Vxc_loc[i];
                Bx_t  =  Bx_loc[i];
                By_t  =  By_loc[i];
                Bz_t  =  Bz_loc[i];

                Vxc_loc[i] = Vxc_t +   Bz_t;
                Bx_loc[i]  = Bx_t  - I*By_t;
                By_loc[i]  = Bx_t  + I*By_t;
                Bz_loc[i]  = Vxc_t -   Bz_t;
            }
        }
        /* Compute elph-matrix elements*/
        for (ND_int i =0 ; i <lattice->kmap->dims[0]; ++i )
        {   
            int ik    = *(kmap + i*2)      ;
            int ksym  = *(kmap + i*2 + 1)  ;
            int ikq   = *(kmap + KplusQidxs[i]*2)      ;
            int kqsym = *(kmap + KplusQidxs[i]*2 + 1)  ;
            ELPH_cmplx * elph_kq_mn = elph_kq + i*elph_kstride_k + iv*elph_kstride_mode ;

            elphLocal(qpt, wfcs, lattice, ikq, ik, kqsym, ksym, true, &dVscf, elph_kq_mn);
        }
    }

    free(KplusQidxs);
    ND_function(uninit,Nd_cmplxS)(&dVscf);
    ND_function(uninit,Nd_cmplxS)(&VlocG_iv);
    free_wfcBox(&wfcr_buffer);
    ND_function(destroy,Nd_cmplxS)(&VlocG);
}

