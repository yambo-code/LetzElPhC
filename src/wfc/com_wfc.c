/*
This file contains functions which gets the wavefunction located 
in another kpool.
*/
#include "wfc.h"

static int MPI_Send_cmplx(const ELPH_cmplx *buf, ND_int count, \
                            int dest, int * tag, MPI_Comm comm);

static int MPI_Send_float(const ELPH_float *buf, ND_int count, \
                            int dest, int * tag, MPI_Comm comm);

static int MPI_Recv_cmplx(const ELPH_cmplx *buf, ND_int count, \
                            int source, int * tag, MPI_Comm comm);

static int MPI_Recv_float(const ELPH_float *buf, ND_int count, \
                            int source, int * tag, MPI_Comm comm);

static void wfc_point2point(struct WFC * wfcs, int ik_send_loc, struct WFC * wfc_out, \
                        int send_rank, int recv_rank, MPI_Comm commQ);




void get_wfc_from_pool(struct WFC * wfcs, ND_int ik, ND_int niBZ_tot,
            MPI_Comm commK, MPI_Comm commQ, struct WFC * wfc_out)
{   
    /*
    This function communicates wfc between kpools.

    ik is global index of k
    */

    int cpus_qpool, qcolor, qrank;

    MPI_Comm_rank(commQ, &qrank);
    MPI_Comm_size(commQ, &cpus_qpool);

    int npw_cpus, krank;
    MPI_Comm_rank(commK, &krank);
    MPI_Comm_size(commK, &npw_cpus);

    int kcolor = qrank/npw_cpus;

    int nkpools = cpus_qpool/npw_cpus ;
    if (nkpools == 1)
    {
        /*
        just copy data to to the out wfc
        */
        struct WFC * wfc_in_tmp = wfcs + ik;

        /*
        Allocate memory for the output wfc. This must be free later
        */
        ND_function(init, Nd_cmplxS) (wfc_out->wfc,  wfc_in_tmp->wfc->rank[0],  wfc_in_tmp->wfc->dims);
        ND_function(init, Nd_floatS) (wfc_out->gvec, wfc_in_tmp->gvec->rank[0], wfc_in_tmp->gvec->dims);
        ND_function(init, Nd_floatS) (wfc_out->Fk,   wfc_in_tmp->Fk->rank[0],   wfc_in_tmp->Fk->dims);

        ND_function(malloc, Nd_cmplxS) (wfc_out->wfc);
        ND_function(malloc, Nd_floatS) (wfc_out->gvec);
        ND_function(malloc, Nd_floatS) (wfc_out->Fk);

        ND_function(copy, Nd_cmplxS) (wfc_in_tmp->wfc, wfc_out->wfc);
        ND_function(copy, Nd_floatS) (wfc_in_tmp->gvec, wfc_out->gvec);
        ND_function(copy, Nd_floatS) (wfc_in_tmp->Fk , wfc_out->Fk);

        wfc_out->npw_total     = wfc_in_tmp->npw_total; 

        wfc_out->npw_loc       = wfc_in_tmp->npw_loc;  

        return ;
    }
    /* get the kcolor that has the required kpoint */
    else    
    {   
        int nks_per_cpu = niBZ_tot/nkpools;
        int nks_rem = niBZ_tot%nkpools;
        /*
        get the kpool where wfc is located
        */
        int in_which_pool = -1;
        int wfc_idx_in_pool = -1;
        for (int ipl = 0 ; ipl < nkpools; ++ipl)
        {
            int kshift = nks_per_cpu*ipl ;
            if (ipl < nks_rem) kshift += ipl ;
            else kshift += nks_rem ;
            int idx_kpool = kshift-ik ;
            if (idx_kpool >= 0)
            {
                in_which_pool = ipl-1;
                wfc_idx_in_pool = -wfc_idx_in_pool;
                break;
            }
            wfc_idx_in_pool = idx_kpool;
        }
        if (in_which_pool < 0 || wfc_idx_in_pool < 0) error_msg("some thing wrong in kpool distribution of wfcs");

        for (int ipwcpu = 0 ; ipwcpu < npw_cpus; ++ipwcpu)
        {
            if (ipwcpu == krank) wfc_point2point(wfcs, wfc_idx_in_pool, wfc_out, \
                                    in_which_pool*npw_cpus+ipwcpu, qrank, commQ);
        }
        
    }
}



static void wfc_point2point(struct WFC * wfcs, int ik_send_loc, struct WFC * wfc_out, \
                        int send_rank, int recv_rank, MPI_Comm commQ)
{   /*
    This function internally allocated large array, that must be freed after this function call
    */
    int qrank;
    MPI_Comm_rank(commQ, &qrank);
    /* send */
    if (qrank == send_rank)
    {   
        struct WFC * wfc_in = wfcs + ik_send_loc;
        int mpi_error;
        /*
        First send basic dimensions of wfc, gvec, Fk arrays
        */
        ND_int rank  = wfc_in->wfc->rank[0];
        ND_int *dims = wfc_in->wfc->dims;
        if (rank > 24 ) error_msg("Max dim supported is 24 \n");

        mpi_error = MPI_Send(&rank, 1, ELPH_MPI_ND_INT, recv_rank, -1, commQ);
        mpi_error = MPI_Send(dims, rank, ELPH_MPI_ND_INT, recv_rank, -2, commQ);

        //gvecs
        rank  = wfc_in->gvec->rank[0];
        dims = wfc_in->gvec->dims;
        if (rank > 24 ) error_msg("Max dim supported is 24 \n");

        mpi_error = MPI_Send(&rank, 1, ELPH_MPI_ND_INT, recv_rank, -1, commQ);
        mpi_error = MPI_Send(dims, rank, ELPH_MPI_ND_INT, recv_rank, -2, commQ);

        //Fk
        rank  = wfc_in->Fk->rank[0];
        dims = wfc_in->Fk->dims;
        if (rank > 24 ) error_msg("Max dim supported is 24 \n");

        mpi_error = MPI_Send(&rank, 1, ELPH_MPI_ND_INT, recv_rank, -1, commQ);
        mpi_error = MPI_Send(dims, rank, ELPH_MPI_ND_INT, recv_rank, -2, commQ);

        //
        ND_int wfc_size = ND_function(size, Nd_cmplxS) (wfc_in->wfc);
        ND_int gvec_size = ND_function(size, Nd_floatS) (wfc_in->gvec);
        ND_int Fk_size = ND_function(size, Nd_floatS) (wfc_in->Fk);
        
        
        mpi_error = MPI_Send(wfc_in->npw_total, 1, ELPH_MPI_ND_INT, recv_rank, 0, commQ);
        mpi_error = MPI_Send(wfc_in->npw_loc, 1, ELPH_MPI_ND_INT, recv_rank, 1, commQ);
        
        int tag = 3 ; // WARNING : this must be same in the recv section too
        mpi_error = MPI_Send_cmplx(wfc_in->wfc->data, wfc_size, recv_rank, &tag, commQ);
        mpi_error = MPI_Send_float(wfc_in->gvec->data, gvec_size, recv_rank, &tag, commQ);
        mpi_error = MPI_Send_float(wfc_in->Fk->data, Fk_size, recv_rank, &tag, commQ);
    }
    else if (qrank == recv_rank)
    {   

        int mpi_error;
        
        ND_int rank;
        ND_int dims[24];
        mpi_error = MPI_Recv(&rank, 1, ELPH_MPI_ND_INT, send_rank, -1, commQ);
        mpi_error = MPI_Recv(dims, rank, ELPH_MPI_ND_INT, send_rank, -2, commQ);
        /*
        Allocate memory for wfc
        */
        ND_function(init, Nd_cmplxS) (wfc_out->wfc,   rank,   dims);
        ND_function(malloc, Nd_cmplxS) (wfc_out->wfc);

        // get rank and dims for gvecs
        mpi_error = MPI_Recv(&rank, 1, ELPH_MPI_ND_INT, send_rank, -1, commQ);
        mpi_error = MPI_Recv(dims, rank, ELPH_MPI_ND_INT, send_rank, -2, commQ);
        /*
        Allocate memory for gvecs
        */
        ND_function(init, Nd_floatS) (wfc_out->gvec,   rank,   dims);
        ND_function(malloc, Nd_floatS) (wfc_out->gvec);
        // get rank and dims for Fk
        mpi_error = MPI_Recv(&rank, 1, ELPH_MPI_ND_INT, send_rank, -1, commQ);
        mpi_error = MPI_Recv(dims, rank, ELPH_MPI_ND_INT, send_rank, -2, commQ);
        /*
        Allocate memory for Fks
        */
        ND_function(init, Nd_floatS) (wfc_out->Fk,   rank,   dims);
        ND_function(malloc, Nd_floatS) (wfc_out->Fk);

        
        ND_int wfc_size = ND_function(size, Nd_cmplxS) (wfc_out->wfc);
        ND_int gvec_size = ND_function(size, Nd_floatS) (wfc_out->gvec);
        ND_int Fk_size = ND_function(size, Nd_floatS) (wfc_out->Fk);

        mpi_error = MPI_Recv(wfc_out->npw_total, 1, ELPH_MPI_ND_INT, send_rank, 0, commQ);
        mpi_error = MPI_Recv(wfc_out->npw_loc, 1, ELPH_MPI_ND_INT, send_rank, 1, commQ);

        
        ND_int int_max = ((ND_int)INT_MAX) -10;
        /* note count can exceed max value of int */
        int tag = 3 ; // WARNING : this must be same in the send section too
        mpi_error = MPI_Recv_cmplx(wfc_out->wfc->data, wfc_size, send_rank, &tag, commQ);
        mpi_error = MPI_Recv_float(wfc_out->gvec->data, gvec_size, send_rank, &tag, commQ);
        mpi_error = MPI_Recv_float(wfc_out->Fk->data, Fk_size, send_rank, &tag, commQ);
    }
}



static int MPI_Send_cmplx(const ELPH_cmplx *buf, ND_int count, \
                            int dest, int * tag, MPI_Comm comm)
{
    /*
    This is just MPI_Send but can send large amounts of data. The purpose is that
    int can overflow
    */
    ND_int max_int_val = ((ND_int)INT_MAX) -10 ;

    int tag_temp = *tag ;
    
    int mpi_error;
    ND_int nsends = count/max_int_val ;
    ND_int rem = count%max_int_val ;

    for (ND_int isends = 0; isends<nsends; ++isends)
    {
        mpi_error = MPI_Send(buf + isends*max_int_val, max_int_val, ELPH_MPI_cmplx, dest, tag_temp, comm);
        ++tag_temp;
    }
    if (rem != 0)
    {
        mpi_error = MPI_Send(buf + nsends*max_int_val, rem, ELPH_MPI_cmplx, dest, tag_temp, comm);
        ++tag_temp;
    }
    *tag = tag_temp;

    return mpi_error;
}



static int MPI_Send_float(const ELPH_float *buf, ND_int count, \
                            int dest, int * tag, MPI_Comm comm)
{
    /*
    This is just MPI_Send but can send large amounts of data. The purpose is that
    int can overflow
    */
    ND_int max_int_val = ((ND_int)INT_MAX) -10 ;

    int tag_temp = *tag ;
    
    int mpi_error;
    ND_int nsends = count/max_int_val ;
    ND_int rem = count%max_int_val ;

    for (ND_int isends = 0; isends<nsends; ++isends)
    {
        mpi_error = MPI_Send(buf + isends*max_int_val, max_int_val, ELPH_MPI_float, dest, tag_temp, comm);
        ++tag_temp;
    }
    if (rem != 0)
    {
        mpi_error = MPI_Send(buf + nsends*max_int_val, rem, ELPH_MPI_float, dest, tag_temp, comm);
        ++tag_temp;
    }
    *tag = tag_temp;

    return mpi_error;
}




static int MPI_Recv_cmplx(const ELPH_cmplx *buf, ND_int count, \
                            int source, int * tag, MPI_Comm comm)
{
    /*
    This is just MPI_Send but can send large amounts of data. The purpose is that
    int can overflow
    */
    ND_int max_int_val = ((ND_int)INT_MAX) -10 ;

    int tag_temp = *tag ;
    
    int mpi_error;
    ND_int nsends = count/max_int_val ;
    ND_int rem = count%max_int_val ;

    for (ND_int isends = 0; isends<nsends; ++isends)
    {
        mpi_error = MPI_Recv(buf + isends*max_int_val, max_int_val, ELPH_MPI_cmplx, source, tag_temp, comm);
        ++tag_temp;
    }
    if (rem != 0)
    {
        mpi_error = MPI_Recv(buf + nsends*max_int_val, rem, ELPH_MPI_cmplx, source, tag_temp, comm);
        ++tag_temp;
    }
    *tag = tag_temp;

    return mpi_error;
}



static int MPI_Recv_float(const ELPH_float *buf, ND_int count, \
                            int source, int * tag, MPI_Comm comm)
{
    /*
    This is just MPI_Send but can send large amounts of data. The purpose is that
    int can overflow
    */
    ND_int max_int_val = ((ND_int)INT_MAX) -10 ;

    int tag_temp = *tag ;
    
    int mpi_error;
    ND_int nsends = count/max_int_val ;
    ND_int rem = count%max_int_val ;

    for (ND_int isends = 0; isends<nsends; ++isends)
    {
        mpi_error = MPI_Recv(buf + isends*max_int_val, max_int_val, ELPH_MPI_float, source, tag_temp, comm);
        ++tag_temp;
    }
    if (rem != 0)
    {
        mpi_error = MPI_Recv(buf + nsends*max_int_val, rem, ELPH_MPI_float, source, tag_temp, comm);
        ++tag_temp;
    }
    *tag = tag_temp;

    return mpi_error;
}


