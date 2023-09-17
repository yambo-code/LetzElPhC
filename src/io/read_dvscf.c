/*
This file contains routines that reads eigen vectors and dvscf
*/
#include "io.h"

void get_FFT_dims(char * file_name, ND_int * nq, ND_int * fft_dims)
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



void read_dvscfq(char * file_name, ND_array(Nd_cmplxS) * eigVec,  \
                ND_array(Nd_cmplxS) *dVscf, ND_int iq, MPI_Comm commK)
{
    
}