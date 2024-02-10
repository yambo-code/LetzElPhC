/* This file contains transpose routines need for ffts*/

#include "fft.h"

void fwd_transpose(struct ELPH_fft_plan * plan)
{
    /* This performs a transpose from 
    (nGxy,Nz_loc) -> (nGxy_loc,Nz) 
    The input and output are stored in plan->nz_buf
    plan->fft_buf is used as scratch space
    */

    /*
    Perform Alltoallv plan->nz_buf -> plan->fft_buf
    */
    int mpi_error;
    int ncpus;
    mpi_error = MPI_Comm_size(plan->comm, &ncpus);

    int * gxy_counts = plan->comm_bufs ;
    int * xy_disp  = gxy_counts + ncpus;
    int * z_counts = gxy_counts + 2*ncpus;
    int * z_disp   = gxy_counts + 3*ncpus;

    mpi_error = MPI_Alltoallv(plan->nz_buf, gxy_counts, xy_disp, \
                            ELPH_MPI_cmplx, plan->fft_data, z_counts, \
                            z_disp, ELPH_MPI_cmplx, plan->comm);
    
    /* Now plan->fft_data has the data. 
    we need to redistribute and store to plan->nz_buf */
    for (ND_int ig = 0 ; ig < plan->nGxyloc; ++ig)
    {   
        ELPH_cmplx * dest_ptr = plan->nz_buf + ig*plan->fft_dims[2];
        for (ND_int i = 0 ; i<ncpus; ++i)
        {   
            ND_int zcount_tmp = z_counts[i]/plan->nGxyloc;
            ND_int z_disp_tmp = z_disp[i]/plan->nGxyloc;
            ELPH_cmplx * src_ptr = plan->fft_data + z_disp[i] + ig*zcount_tmp;
            memcpy(dest_ptr+z_disp_tmp, src_ptr, sizeof(ELPH_cmplx)*zcount_tmp);
        }
    }

}




void bwd_transpose(struct ELPH_fft_plan * plan)
{
    /* This performs a transpose from 
    (nGxy_loc,Nz)  -> (nGxy,Nz_loc)
    The input and output are stored in plan->nz_buf
    plan->fft_buf is used as scratch space
    */

    /*
    Perform Alltoallv plan->nz_buf -> plan->fft_buf
    */
    int mpi_error;
    int ncpus;
    mpi_error = MPI_Comm_size(plan->comm, &ncpus);

    int * gxy_counts = plan->comm_bufs ;
    int * xy_disp  = gxy_counts + ncpus;
    int * z_counts = gxy_counts + 2*ncpus;
    int * z_disp   = gxy_counts + 3*ncpus;
    
    /* 
    we need to redistribute and store to plan->fft_data */
    for (ND_int ig = 0 ; ig < plan->nGxyloc; ++ig)
    {   
        ELPH_cmplx * src_ptr = plan->nz_buf + ig*plan->fft_dims[2];
        for (ND_int i = 0 ; i<ncpus; ++i)
        {   
            ND_int zcount_tmp = z_counts[i]/plan->nGxyloc;
            ND_int z_disp_tmp = z_disp[i]/plan->nGxyloc;
            ELPH_cmplx * dest_ptr = plan->fft_data + z_disp[i] + ig*zcount_tmp;
            memcpy(dest_ptr, src_ptr+z_disp_tmp, sizeof(ELPH_cmplx)*zcount_tmp);
        }
    }

    // Now plan->fft_data has the data. scatter back to plan->nz_buf
    mpi_error = MPI_Alltoallv(plan->fft_data, z_counts, z_disp, \
                            ELPH_MPI_cmplx, plan->nz_buf, gxy_counts, \
                            xy_disp, ELPH_MPI_cmplx, plan->comm);

}


