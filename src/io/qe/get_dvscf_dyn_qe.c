/*
This function gets the dvscf and eig vectors for q point from qe
*/
#include "qe_io.h"

void get_dvscf_dyn_qe(const char * ph_save_dir, struct Lattice * lattice, \
        ND_int iq_BZ, ELPH_cmplx * eig, ELPH_cmplx * dvscf, ELPH_float * omega_ph, \
        const struct ELPH_MPI_Comms * Comm)
{   

    // ph_save_dir must be available on all cpus else U.B
    // if dvscf == NULL, then only eig vecs are returned.
    ND_int nmodes = lattice->atomic_pos->dims[0]*3;

    char * tmp_char_buf = malloc(1024+strlen(ph_save_dir));
    ELPH_cmplx * pat_vecs = NULL;
    if (dvscf != NULL) pat_vecs = malloc(sizeof(ELPH_cmplx)*nmodes*nmodes);
    
    if (Comm->commQ_rank ==0)
    {   
        if (dvscf != NULL)
        {
            sprintf(tmp_char_buf,"%s/patterns.%d.xml",ph_save_dir, (int)(iq_BZ+1));
            read_pattern_qe(tmp_char_buf, lattice, pat_vecs);
        }
        ELPH_float qpts[3];
        sprintf(tmp_char_buf,"%s/dyn%d",ph_save_dir, (int)(iq_BZ+1));
        ND_int ndyn_read = read_dyn_qe(tmp_char_buf, lattice, qpts, omega_ph, eig);
        if(ndyn_read != 1) error_msg("Wrong number of dynmats read from dynamical matrix file");
    }

    // Bcast variables
    
    MPI_Bcast(eig,nmodes*nmodes,ELPH_MPI_cmplx,0,Comm->commQ);
    MPI_Bcast(omega_ph,nmodes,ELPH_MPI_float,0,Comm->commQ);

    if (dvscf != NULL)
    {   
        MPI_Bcast(pat_vecs,nmodes*nmodes,ELPH_MPI_cmplx,0,Comm->commQ);

        sprintf(tmp_char_buf,"%s/dvscf%d",ph_save_dir, (int)(iq_BZ+1));
        if (Comm->commRq_rank == 0)
            read_dvscf_qe(tmp_char_buf, lattice, eig, pat_vecs, dvscf, Comm->commK);
        // broad cast dvscf
        ND_int dvscf_len  = nmodes*lattice->nmag*lattice->fft_dims[0]*lattice->fft_dims[1]*lattice->nfftz_loc;
        ND_int max_int_val = ((ND_int)INT_MAX) -10 ;

        // check for buffer overflow
        if (dvscf_len > max_int_val) error_msg("Buffer Overflow in dvscf_len. Contact developer");
        MPI_Bcast(dvscf, dvscf_len, ELPH_MPI_cmplx, 0, Comm->commRq);
        free(pat_vecs);
    }

    free(tmp_char_buf);
    

}


