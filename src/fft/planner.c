#include "fft.h"
#include "string.h"
#include "../common/parallel.h"
/*
Create plan creation function for FFTs
*/

void wfc_plan(struct ELPH_fft_plan * plan, const ND_int ngvecs_loc, const ND_int nzloc, \
            const ND_int nGxyloc, const int * gvecs, const ND_int * fft_dims, \
            unsigned fft_flags, MPI_Comm comm)
{   
    int ncpus, mpi_error,my_rank; 
    mpi_error = MPI_Comm_size(comm, &ncpus);
    mpi_error = MPI_Comm_rank(comm, &my_rank);

    plan->comm_bufs = malloc(sizeof(int)*4*ncpus);
    /*
    gxy_counts, xy_disp, z_counts, z_disp
    */

    plan->nzloc = get_mpi_local_size_idx(fft_dims[2], NULL,  comm);
    if (plan->nzloc != nzloc) error_msg("Wrong z local dimensions to planner.");


    ND_array(Nd_cmplxS) fft_data[1] ;
    ND_function(init, Nd_cmplxS) (fft_data, 3, nd_idx{fft_dims[0],fft_dims[1],plan->nzloc});


    ND_int size_fft_data = fft_dims[0]*fft_dims[1]*plan->nzloc ;
    if (size_fft_data < (fft_dims[2]*nGxyloc) ) size_fft_data = fft_dims[2]*nGxyloc ; 

    plan->align_len = alignment_len();
    size_fft_data += plan->align_len; // we add alignment_len to make alignment_len plans for x and y

    // alloc memory for fft_data and nz_buf;
    plan->fft_data = fftw_fun(malloc)(size_fft_data*sizeof(ELPH_cmplx));
    plan->nz_buf   = fftw_fun(malloc)(size_fft_data*sizeof(ELPH_cmplx)); 

    if (plan->fft_data == NULL) error_msg("fft buffer allocation failed");
    if (plan->nz_buf == NULL)   error_msg("nz  buffer allocation failed");
    
    memset(plan->fft_data,0,size_fft_data*sizeof(ELPH_cmplx)); // 0 the buffer
    memset(plan->nz_buf,  0,size_fft_data*sizeof(ELPH_cmplx)); // 0 the buffer

    plan->comm        = comm;
    plan->gvecs       = gvecs;
    plan->nGxyloc     = nGxyloc;
    plan->ngvecs_loc  = ngvecs_loc;

    // fft_dims Nx,Ny,Nz
    memcpy(plan->fft_dims,fft_dims,sizeof(ND_int)*3);

    mpi_error = MPI_Allreduce(&nGxyloc, &(plan->nGxy), 1, ELPH_MPI_ND_INT, MPI_SUM, comm);

    if (nGxyloc != get_mpi_local_size_idx(plan->nGxy, NULL,comm)) error_msg("Wrong xy local dimensions to planner.");
    
    int * Gxy_loc   = malloc(sizeof(int)*2*nGxyloc);       // local gxy (nGxyloc,2)
    plan->Gxy_total = malloc(sizeof(int)*2*plan->nGxy);

    int * gxy_counts = plan->comm_bufs ;
    int * xy_disp  = gxy_counts + ncpus;
    int * z_counts = gxy_counts + 2*ncpus;
    int * z_disp   = gxy_counts + 3*ncpus;

    int xycount = 0;
    // initialize with gvecs which are out of the box
    int Gx_prev = fft_dims[0]+10;
    int Gy_prev = fft_dims[1]+10;

    plan->ngxy_z = malloc(sizeof(int)*nGxyloc);
    // zero the buffer
    for (ND_int i = 0 ; i < nGxyloc; ++i) plan->ngxy_z[i] = 0;

    int Gxmin = gvecs[0] ; int Gxmax = gvecs[0]; // min and max values of Gx

    for (ND_int ig = 0; ig < plan->ngvecs_loc; ++ig)
    {   
        int Gx_temp = gvecs[3*ig];
        if (Gx_temp > 0) Gx_temp = get_miller_idx(Gx_temp, fft_dims[0]);
        if (Gx_temp < Gxmin) Gxmin = Gx_temp;
        if (Gx_temp > Gxmax) Gxmax = Gx_temp;

        if (gvecs[3*ig+1] != Gy_prev || gvecs[3*ig] != Gx_prev)
        {   
            Gx_prev = gvecs[3*ig];
            Gy_prev = gvecs[3*ig+1];

            Gxy_loc[2*xycount]   = Gx_prev;
            Gxy_loc[2*xycount+1] = Gy_prev;

            if (Gx_prev < 0)   Gxy_loc[2*xycount]   += fft_dims[0];
            if (Gy_prev < 0)   Gxy_loc[2*xycount+1] += fft_dims[1];

            ++xycount;
        }
        
        plan->ngxy_z[xycount-1] += 1;
    }

    if (xycount != nGxyloc) error_msg("Wrong number of xy components in each cpu");
    // get global Gmin anf Gmax
    mpi_error = MPI_Allreduce(MPI_IN_PLACE, &Gxmin, 1, MPI_INT, MPI_MIN, comm);
    mpi_error = MPI_Allreduce(MPI_IN_PLACE, &Gxmax, 1, MPI_INT, MPI_MAX, comm);
    // so we have to perform ffts from [0,Gxmax] and [Nx-Gxmin,Gxmax]


    mpi_error = MPI_Allgather(&xycount, 1, MPI_INT, gxy_counts, 1, MPI_INT, comm);

    int nzloc_temp = nzloc;
    mpi_error = MPI_Allgather(&nzloc_temp, 1, MPI_INT, z_counts, 1, MPI_INT, comm);

    int dispdisp_temp = 0;
    for (int i = 0 ; i <ncpus; ++i)
    {   
        gxy_counts[i] *= 2 ;
        xy_disp[i] = dispdisp_temp;
        dispdisp_temp += gxy_counts[i];
    }

    mpi_error = MPI_Allgatherv(Gxy_loc, 2*xycount, MPI_INT, plan->Gxy_total, gxy_counts, \
                                xy_disp, MPI_INT, comm);

    free(Gxy_loc);

    /* fill the comm buffer */
    dispdisp_temp = 0;
    for (int i = 0 ; i <ncpus; ++i)
    {
        gxy_counts[i] = nzloc*(gxy_counts[i]/2);
        xy_disp[i]    = nzloc*(xy_disp[i]/2);

        z_counts[i] *= xycount;
        z_disp[i] = dispdisp_temp;
        dispdisp_temp += z_counts[i];
    }

    /* Now time to create plans */
    // i) create plan along entire X direction
    // ii) create two plans along y for -Gx_min < Gx < Gx_max
    // iii) create plan for along z only for set of (Gx,Gy) pairs
    
    // create plan buffers
    plan->fplan_x  = malloc(sizeof(ND_function(FFT_plan, Nd_cmplxS))* 6 * plan->align_len);  // (naligment plans) for x 
    plan->fplan_y  = plan->fplan_x + plan->align_len;  // (naligment plans) for y
    plan->fplan_y1 = plan->fplan_x + 2*plan->align_len; // (naligment plans) for y

    /* backward plans G->r */
    plan->bplan_x  = plan->fplan_x + 3*plan->align_len;  // (naligment plans) for x 
    plan->bplan_y  = plan->fplan_x + 4*plan->align_len;  // (naligment plans) for y
    plan->bplan_y1 = plan->fplan_x + 5*plan->align_len; // (naligment plans) for y
    
    
    if (Gxmin >= 0 || Gxmax <= 0) 
    {   
        /* In this case something is fishy, so do a full fft along x,y */
        if (fft_dims[0]%2 == 0)
        {
            Gxmin = -fft_dims[0]/2;
            Gxmax =  (fft_dims[0]/2)-1;
        }
        else
        {
            Gxmin = -(fft_dims[0]-1)/2;
            Gxmax =  (fft_dims[0]-1)/2; 
        }
    }
    
    Gxmin += fft_dims[0];    
    plan->Gxmin = Gxmin;   
    plan->Gxmax = Gxmax;

    // DEbug// comment me later
    if (Gxmin <0 || Gxmax <0) error_msg("FFT fail along xy");

    // i) create forward plan and bwd plan along X
    for (ND_int i = 0; i<plan->align_len; ++i)
    {
        fft_data->data = plan->fft_data+i;
        const ND_int in_idx[2] = {0,1};
        
        ND_int ia = fftw_fun(alignment_of)((void *)fft_data->data);
        ia /= sizeof(ELPH_cmplx);

        ND_function(fft_planner, Nd_cmplxS) (fft_data, fft_data, 1, \
                    in_idx, -1, NULL, fft_flags, plan->fplan_x+ia);
        if (plan->fplan_x[ia] == NULL) error_msg("X forward plan failed");

        // backward plan G->r
        ND_function(fft_planner, Nd_cmplxS) (fft_data, fft_data, 1, \
                    in_idx, +1, NULL, fft_flags, plan->bplan_x+ia);
        if (plan->bplan_x[ia] == NULL) error_msg("X backward plan failed");
    }

    // uninitiate the fft_data
    ND_function(uninit, Nd_cmplxS) (fft_data);
    ND_function(init, Nd_cmplxS) (fft_data, 3, nd_idx{(Gxmax+1),fft_dims[1],plan->nzloc});

    
    // ii) create two Y ffts plans along Y :  N-Gmin < Gx < N-1 and 

    //  first we create Y plans for 0 <= Gx <= Gmax
    for (ND_int i = 0; i<plan->align_len; ++i)
    {
        fft_data->data = plan->fft_data+i;
        const ND_int in_idx = 1;
        
        ND_int ia = fftw_fun(alignment_of)((void *)fft_data->data);
        ia /= sizeof(ELPH_cmplx);

        ND_function(fft_planner, Nd_cmplxS) (fft_data, fft_data, 1, \
                    &in_idx, -1, NULL, fft_flags, plan->fplan_y+ia);
        if (plan->fplan_y[ia] == NULL) error_msg("Y forward plan failed");

        // backward plan G->r
        ND_function(fft_planner, Nd_cmplxS) (fft_data, fft_data, 1, \
                    &in_idx, +1, NULL, fft_flags, plan->bplan_y+ia);
        if (plan->bplan_y[ia] == NULL) error_msg("Y backward plan failed");
    }
    
    ND_function(uninit, Nd_cmplxS) (fft_data);
    ND_function(init, Nd_cmplxS) (fft_data, 3, nd_idx{(fft_dims[0]-Gxmin),fft_dims[1],plan->nzloc});


    // Then we create Y plans for N-Gmin < Gx < N-1
    for (ND_int i = 0; i<plan->align_len; ++i)
    {
        fft_data->data = plan->fft_data+i;
        const ND_int in_idx = 1;
        
        ND_int ia = fftw_fun(alignment_of)((void *)fft_data->data);
        ia /= sizeof(ELPH_cmplx);

        ND_function(fft_planner, Nd_cmplxS) (fft_data, fft_data, 1, \
                    &in_idx, -1, NULL, fft_flags, plan->fplan_y1+ia);
        if (plan->fplan_y1[ia] == NULL) error_msg("Y1 forward plan failed");

        // backward plan G->r
        ND_function(fft_planner, Nd_cmplxS) (fft_data, fft_data, 1, \
                    &in_idx, +1, NULL, fft_flags, plan->bplan_y1+ia);
        if (plan->bplan_y1[ia] == NULL) error_msg("Y1 backward plan failed");
    }
    ND_function(uninit, Nd_cmplxS) (fft_data);

    // iii) create a single z plan. 
    // forward plan-> 
    plan->fplan_z = fftw_fun(plan_many_dft)(1, (int[1]){fft_dims[2]}, nGxyloc, plan->nz_buf, \
                        NULL, 1, fft_dims[2], plan->nz_buf, NULL,1, fft_dims[2],-1, fft_flags);
    if (plan->fplan_z == NULL) error_msg("Z forward plan failed");
    
    // backward plan 
    plan->bplan_z = fftw_fun(plan_many_dft)(1, (int[1]){fft_dims[2]}, nGxyloc, plan->nz_buf, \
                        NULL, 1, fft_dims[2], plan->nz_buf, NULL,1, fft_dims[2],+1, fft_flags);
    if (plan->bplan_z == NULL) error_msg("Z backward plan failed");

}


void wfc_destroy_plan(struct ELPH_fft_plan * plan)
{   
    // destroy x,y,z buffer
    fftw_fun(free)(plan->fft_data);
    fftw_fun(free)(plan->nz_buf);
    // destroy plans

    for (ND_int i = 0; i<plan->align_len; ++i)
    {
        ND_function(fft_destroy_plan, Nd_cmplxS) (plan->fplan_x[i]);
        ND_function(fft_destroy_plan, Nd_cmplxS) (plan->bplan_x[i]);

        ND_function(fft_destroy_plan, Nd_cmplxS) (plan->fplan_y[i]);
        ND_function(fft_destroy_plan, Nd_cmplxS) (plan->bplan_y[i]);
        ND_function(fft_destroy_plan, Nd_cmplxS) (plan->fplan_y1[i]);
        ND_function(fft_destroy_plan, Nd_cmplxS) (plan->bplan_y1[i]);
        
    }

    ND_function(fft_destroy_plan, Nd_cmplxS) (plan->fplan_z);
    ND_function(fft_destroy_plan, Nd_cmplxS) (plan->bplan_z);

    // free remaining buffers
    free(plan->fplan_x);
    free(plan->comm_bufs);
    free(plan->Gxy_total);
    free(plan->ngxy_z);

}




