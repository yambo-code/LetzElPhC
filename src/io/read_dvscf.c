/*
This file contains routines that reads eigen vectors and dvscf
*/
#include "io.h"

void get_FFT_dims(const char * file_name, ND_int * nq, ND_int * fft_dims)
{   

    int my_rank, comm_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);

    int Nfft[3];
    int nqpts;
    if (my_rank == 0)
    {
        int ncid, varid, retval, nc_rank; // variables 

        if ((retval = nc_open(file_name, NC_NOWRITE, &ncid))) ERR(retval); // open file

        /* read fft dims */
        if ((retval = nc_inq_varid(ncid, "fft_dims", &varid))) ERR(retval); // get the varible id of the file
        if ((retval = nc_get_var(ncid, varid, Nfft))) ERR(retval); //get data in floats
        /* read number of q qpoints */
        if ((retval = nc_inq_varid(ncid, "qpts", &varid))) ERR(retval); // get the varible id of the file
        if ((retval = nc_get_var(ncid, varid, &nqpts))) ERR(retval); //get data in floats

        /*close the file*/
        if ((retval = nc_close(ncid))) ERR(retval); // close the file
    }

    MPI_Bcast(Nfft,   3, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&nqpts, 1, MPI_INT, 0, MPI_COMM_WORLD);

    *nq         = nqpts;
    fft_dims[0] = Nfft[0];
    fft_dims[1] = Nfft[1];
    fft_dims[2] = Nfft[2];
    

}



void read_dvscfq(const char * file_name, ND_array(Nd_cmplxS) * eigVec,  \
                struct Lattice * lattice, ND_array(Nd_cmplxS) *dVscf, \
                ND_int iq, MPI_Comm commK)
{
    /*
    // WARNING :  memory is allocated internally, must be free outside the function
    */
    int retval, dvscfid,var_id;

    int my_rank;
    MPI_Comm_rank(commK, &my_rank);
    
    size_t startp[5] = {iq,0,0,lattice->nfft_shift_loc,0} ;
    size_t countp[5];
    int dim_ids[5];
    // ('nqpts','nmodes','nmag','nfft','re_im')
    if ((retval = nc_open_par(file_name, NC_NOWRITE, commK, MPI_INFO_NULL, &dvscfid))) ERR(retval);
    if ((retval = nc_inq_varid(dvscfid, "dVscfs", &var_id))) ERR(retval); // get the id of the req variable
    if ((retval = nc_inq_var(dvscfid, var_id, NULL, NULL, NULL, dim_ids, NULL ))) ERR(retval); // get dims
    //
    for (ND_int i = 0; i < 5; ++i)
    {
        if ((retval = nc_inq_dimlen(dvscfid, dim_ids[i], countp + i))) ERR(retval);
    }
    countp[0] = 1; countp[3] = lattice->nffts_loc;

    ND_function(init,Nd_cmplxS)(dVscf, 3, nd_idx{countp[1],countp[2],lattice->nffts_loc});
    ND_function(malloc,Nd_cmplxS)(dVscf);

    if ((retval = nc_var_par_access(dvscfid, var_id, NC_COLLECTIVE))) ERR(retval); // NC_COLLECTIVE or NC_INDEPENDENT

    if ((retval = nc_get_vara(dvscfid, var_id, startp, countp, dVscf->data))) ERR(retval); //get data in floats

    /*
    read the phonon modes
    */
    if(my_rank == 0)
    {
        // ('nqpts', 'nmodes', 'natom','xcart','re_im' ))
        size_t startp_ph[5] = {iq,0,0,0,0};
        if ((retval = nc_inq_varid(dvscfid, "ph_pol_vec", &var_id))) ERR(retval); // get the id of the req variable
        if ((retval = nc_inq_var(dvscfid, var_id, NULL, NULL, NULL, dim_ids, NULL ))) ERR(retval); // get dims
        for (ND_int i = 0; i < 5; ++i)
        {
            if ((retval = nc_inq_dimlen(dvscfid, dim_ids[i], countp + i))) ERR(retval);
        }
        countp[0] = 1;
        if ((retval = nc_var_par_access(dvscfid, var_id, NC_INDEPENDENT))) ERR(retval);

        ND_function(init,Nd_cmplxS)(eigVec, 3, nd_idx{countp[1],countp[2],countp[3]});
        ND_function(malloc,Nd_cmplxS)(eigVec);

        if ((retval = nc_get_vara(dvscfid, var_id, startp_ph, countp, eigVec->data))) ERR(retval); //get data in floats

    }
    
    Bcast_ND_arrayCmplx(eigVec, true, 0, commK);

    NC_close_file(dvscfid);

}

