#pragma once

#include <mpi.h>

#include "common/dtypes.h"
#include "elphC.h"

/* This struct contains all the elph calculation input file details */
struct elph_usr_input
{
    // system varibles
    int nkpool;         // k point parallelization
    int nqpool;         // q point parallelization
    int start_bnd;      // starting band
    int end_bnd;        // last band
    char* save_dir;     // save dir
    char* ph_save_dir;  // ph_save directory
    char* kernel_str;   // level of screening to include
    bool kminusq;       // true if convention is "yambo" else false
};

/* This struct contains all the interpolation calculation input file details */
struct interpolation_usr_input
{
    char* ph_save_dir;                // ph_save dir
    char* ph_save_interpolation_dir;  // directory for interpolation files
    bool interpolate_dvscf;           // interpolate dvscf
    bool asr;                         // apply acoustic sum rule
    char asr_kind[32];                // type of asr to be applied.
    bool loto;                        // if true apply LO-TO splitting
    ELPH_float loto_dir[3];           // LO-TO splitting direction
    ND_int qgrid_fine[3];             // fine grid interpolation.
    bool write_dVbare;  // If True, will dump change local part of dVbare
    ELPH_float
        eta_induced;    // ewald parameter for screened long range screened
                        // interactions (dipoles, quadrupoles, phonon dynmats)
    ELPH_float eta_ph;  // ewald summation parameter for interpolation of
                        // dynamical matrices.
    char* qlist_file;   // in case you want to interpolate over a given list,
                        // file name
};

// electron-phonon input file parser
void init_elph_usr_input(struct elph_usr_input** input);
void free_elph_usr_input(struct elph_usr_input* input);
void read_elph_input_file(const char* input_file,
                          struct elph_usr_input** input_data,
                          MPI_Comm MPI_world_comm);

// interpolation input io
void init_interpolation_usr_input(struct interpolation_usr_input** input);
void free_interpolation_usr_input(struct interpolation_usr_input* input);
void read_interpolation_input_file(const char* input_file,
                                   struct interpolation_usr_input** input_data,
                                   MPI_Comm MPI_world_comm);
