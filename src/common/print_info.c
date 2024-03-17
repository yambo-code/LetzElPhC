#include "print_info.h"
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdint.h>
#include <mpi.h>
#include <string.h>
#include <stdarg.h>
#include <stdlib.h>

void print_info_msg(int mpi_rank, const char* fmt, ...)
{
    if (mpi_rank)
    {
        return;
    }
    va_list args;
    va_start(args, fmt);
    vfprintf(stdout, fmt, args);
    va_end(args);
    fprintf(stdout, "\n");
    fflush(stdout);
}

void print_input_info(const char* save_dir, const char* ph_save_dir,
                      const char* kernel, const bool kminusq, const enum ELPH_dft_code dft_code,
                      const struct ELPH_MPI_Comms* Comm)
{
    print_info_msg(Comm->commW_rank, "");
    print_info_msg(Comm->commW_rank, "============== Basic info ==============");
    char convention[16];
    if (kminusq)
    {
        snprintf(convention, sizeof(convention), "%s", "k-q -> k");
    }
    else
    {
        snprintf(convention, sizeof(convention), "%s", "k -> k+q");
    }

    if (dft_code == DFT_CODE_QE)
    {
        print_info_msg(Comm->commW_rank, "DFT Code                : %s", "Quantum ESPRESSO");
    }
    else
    {
        print_info_msg(Comm->commW_rank, "DFT Code                : %s", "Unknown");
    }
    print_info_msg(Comm->commW_rank, "SAVE Dir                : %s", save_dir);
    print_info_msg(Comm->commW_rank, "ph_save Dir             : %s", ph_save_dir);
    print_info_msg(Comm->commW_rank, "Screening Kernel        : %s", kernel);
    print_info_msg(Comm->commW_rank, "ELPH Output convention  : %s", convention);

    print_info_msg(Comm->commW_rank, "");
    print_info_msg(Comm->commW_rank, "============== Parallelization info ==============");
    print_info_msg(Comm->commW_rank, "Total cpus              : %d", Comm->commW_size);
    print_info_msg(Comm->commW_rank, "Number of q-pools       : %d", Comm->nqpools);
    print_info_msg(Comm->commW_rank, "Number of k-pools       : %d", Comm->nkpools);
    print_info_msg(Comm->commW_rank, "Cpus per k-pool         : %d", Comm->commK_size);
#if defined(ELPH_OMP_PARALLEL_BUILD)
    char* omp_n_threads = getenv("OMP_NUM_THREADS");
    if (omp_n_threads)
    {
        print_info_msg(Comm->commW_rank, "Openmp Threads    = %s", omp_n_threads);
    }
    else
    {
        print_info_msg(Comm->commW_rank, "Number of OMP Threads : 1");
    }
#endif
}

void print_lattice_info(const struct ELPH_MPI_Comms* Comm, const struct Lattice* lattice)
{
    print_info_msg(Comm->commW_rank, "");
    print_info_msg(Comm->commW_rank, "============== Lattice info ==============");
    print_info_msg(Comm->commW_rank, "Number of atoms in Unit Cell  : %d", (int)lattice->natom);
    print_info_msg(Comm->commW_rank, "Number of spin components     : %d", (int)lattice->nspin);
    print_info_msg(Comm->commW_rank, "Number of spinor components   : %d", (int)lattice->nspinor);

    print_info_msg(Comm->commW_rank, "Dimension of system           : %cD", lattice->dimension);
    //
    if (lattice->is_soc_present)
    {
        print_info_msg(Comm->commW_rank, "Spin orbit coupling           : True");
    }
    else
    {
        print_info_msg(Comm->commW_rank, "Spin orbit coupling           : False");
    }
    //
    if (lattice->nmag == 1)
    {
        print_info_msg(Comm->commW_rank, "Magnetic calculation type     : Non-magnetic");
    }
    else if (lattice->nmag == 2)
    {
        print_info_msg(Comm->commW_rank, "Magnetic calculation type     : LSDA (collinear spin)");
    }
    else
    {
        print_info_msg(Comm->commW_rank, "Spin orbit coupling           : magnetic non-collinear");
    }
    //
    print_info_msg(Comm->commW_rank, "Number of kpoints in iBZ      : %d", (int)lattice->nkpts_iBZ);
    print_info_msg(Comm->commW_rank, "Number of kpoints in full BZ  : %d", (int)lattice->nkpts_BZ);
    print_info_msg(Comm->commW_rank, "Number of symmetries          : %d", (int)lattice->nsym);
    //
    if (lattice->timerev == 1)
    {
        print_info_msg(Comm->commW_rank, "Time reversal for electrons   : True");
    }
    else
    {
        print_info_msg(Comm->commW_rank, "Time reversal for electrons   : False");
    }
    //
    print_info_msg(Comm->commW_rank, "FFT Grid                      : %d  %d  %d",
                   (int)lattice->fft_dims[0], (int)lattice->fft_dims[1], (int)lattice->fft_dims[2]);

    print_info_msg(Comm->commW_rank, "Bands range for elph calc     : %d - %d",
                   (int)lattice->start_band, (int)lattice->end_band);
}

void print_phonon_info(const struct ELPH_MPI_Comms* Comm, const struct Phonon* phonon)
{
    print_info_msg(Comm->commW_rank, "");
    print_info_msg(Comm->commW_rank, "============== Phonon info ==============");

    print_info_msg(Comm->commW_rank, "Number of qpoints in iBZ      : %d", (int)phonon->nq_iBZ);
    print_info_msg(Comm->commW_rank, "Number of qpoints in full BZ  : %d", (int)phonon->nq_BZ);
    print_info_msg(Comm->commW_rank, "Number of phonon symmetries   : %d", (int)phonon->nph_sym);

    int time_rev_present = 0;
    for (int isym = 0; isym < phonon->nph_sym; ++isym)
    {
        if (phonon->ph_syms[isym].time_rev)
        {
            time_rev_present = 1;
            break;
        }
    }

    if (time_rev_present == 1)
    {
        print_info_msg(Comm->commW_rank, "Time reversal for phonons     : True");
    }
    else
    {
        print_info_msg(Comm->commW_rank, "Time reversal for phonons     : False");
    }
}