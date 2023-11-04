#include "dvloc.h"

/* create a interpolation table on a coarse grid for short local potential in Gspace*/
void vlocg_table(const struct Lattice * lattice, const struct Pseudo * pseudo, 
                ND_int * npts_co, ELPH_float ** g_co, ELPH_float ** vlocg, \
                ELPH_float ** vploc_co, MPI_Comm commK)
{   
    /*
    g_co, vlocg, vploc_co are allocate internally and must be freed outside of this function
    */
    int ncpus;
    int mpi_error = MPI_Comm_size(commK, &ncpus);

    const ND_array(Nd_floatS)* Vloc_atomic  = pseudo->Vloc_atomic;
    const ND_array(Nd_floatS)* r_grid       = pseudo->r_grid;
    const ND_array(Nd_floatS)* rab_grid     = pseudo->rab_grid;
    const ELPH_float * Zval                 = pseudo->Zval;
    const ND_int ntype                      = pseudo->ntype;
    const ND_int ngrid_max                  = pseudo->ngrid_max;
    ELPH_float volume                       = fabs(det3x3(lattice->alat_vec->data));

    ELPH_float gmax, gmin;
    ND_int npts;
    
    // create a small scope
    {   
    ELPH_float blat[9]; 
    reciprocal_vecs(lattice->alat_vec->data,blat); // 2*pi is included 
    
    ELPH_float bi[3], bmin, bmax, Nmax;
    
    bi[0] =  sqrt((blat[0]*blat[0]) + (blat[3]*blat[3]) + (blat[6]*blat[6]));
    bi[1] =  sqrt((blat[1]*blat[1]) + (blat[4]*blat[4]) + (blat[7]*blat[7]));
    bi[2] =  sqrt((blat[2]*blat[2]) + (blat[5]*blat[5]) + (blat[8]*blat[8]));
    bmin = bi[0]; bmax = bi[0];  Nmax = lattice->fft_dims[0];
    
    for (int i = 0 ; i < 3 ; ++i)
    {
        if (bi[i]<bmin) bmin = bi[i];
        if (bi[i]>bmax) bmax = bi[i];
        if (lattice->fft_dims[i]>Nmax) Nmax = lattice->fft_dims[i];
    }

    /*
    We try to do this only once for any q, so we should get the max possible value of |G+q|
    |G+q| <= |G| + |q| 
    if all three values of q (in crystal coordinates) are restricted to [-1,1]
    then |q| <= 3*b_max
    G_max = |N1/2*b1 + N2/2*b2 + N3/2*b3| < N_max/2 * 3 * b_max
    |G|<= 1.5*N_max*b_max
    => |G+q| <= b_max*(1.5*N_max+3)
    So we choose Interpolation range to be [0,(1.5*N_max+3)*b_max]
    // coarse spacing = b_min, npts ~ (1.5*N_max + 3)*b_max/b_min. 
    */

    gmax = (1.7*Nmax +3)*bmax; // we choose 1.7 instead of 1.5
    npts = ceil(gmax/bmin);
    gmin = 0.0;
    }

    
    ELPH_float * xins = malloc(sizeof(ELPH_float)*npts);
    ELPH_float * yins = malloc(sizeof(ELPH_float)*npts*ntype);
    *vploc_co = malloc(sizeof(ELPH_float)*npts*ntype);

    *npts_co = npts;
    *g_co    = xins;
    *vlocg   = yins;

    ELPH_float diff = (gmax-gmin)/(npts-1.0) ;
    for (ND_int i=0; i<npts; ++i) xins[i] = gmin + i*diff ;

    ND_int n_shift;
    int npts_loc = get_mpi_local_size_idx(npts, &n_shift, commK);

    ELPH_float * work_array = malloc(sizeof(ELPH_float)*(ngrid_max+npts_loc));
    ELPH_float * vlocg_cpu  = work_array + ngrid_max;

    int * counts_recv = malloc(sizeof(int)*2*ncpus);
    int * displacements = counts_recv + ncpus;

    mpi_error = MPI_Allgather(&npts_loc, 1, MPI_INT, counts_recv, 1, MPI_INT, commK);
    int n_shift_temp = n_shift;
    mpi_error = MPI_Allgather(&n_shift_temp, 1, MPI_INT, displacements, 1, MPI_INT, commK);
    
    for (ND_int itype = 0; itype <ntype; ++itype )
    {   
        ELPH_float * restrict vlocg_atom = yins + itype*npts;
        ELPH_float * restrict dyy = *vploc_co + itype*npts;
        // loop not thread safe
        for (ND_int i = 0 ; i < npts_loc; ++i)
        {
            vlocg_cpu[i] = Vloc_Gspace(work_array, lattice->dimension, \
                                    xins[i+n_shift], Vloc_atomic+itype, \
                                    r_grid+itype, rab_grid+itype, Zval[itype], 1.0, volume);
        }
        mpi_error = MPI_Allgatherv(vlocg_cpu, npts_loc, ELPH_MPI_float, vlocg_atom, \
                                counts_recv, displacements, ELPH_MPI_float, commK);
        
        prepare_spline(npts, xins, vlocg_atom, dyy);
    }

    free(counts_recv);
    free(work_array);

}


