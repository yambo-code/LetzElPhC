#include "../common/dtypes.h"
#include "../common/error.h"
#include "../elphC.h"
#include "io.h"
#include <netcdf.h>
#include <netcdf_par.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// write basic data related to the dft and ph
//
void write_basic_data(const int ncid, struct Lattice* lattice,
                      struct Phonon* phonon, const char* kernel_str, const char* convention_str)
{
    // only single cpu must call this which implies that the file must be opened
    // by only single cpu
    //
    //  writes kpoints in (crystal units)
    //  qpoints in (crystal units)
    //  kernel type
    //  convention
    //  start and end band
    //  nsym, symmetries matrices for phonons
    int nc_err;
    int varid;

    // kpoints
    def_ncVar(ncid, &varid, 2, ELPH_NC4_IO_FLOAT,
              (ND_int[]) { lattice->nkpts_BZ, 3 }, "kpoints", (char*[]) { "nk", "pol" },
              NULL);
    if ((nc_err = nc_put_var(ncid, varid, lattice->kpt_fullBZ_crys)))
    {
        ERR(nc_err);
    }

    // qpoints
    def_ncVar(ncid, &varid, 2, ELPH_NC4_IO_FLOAT,
              (ND_int[]) { phonon->nq_BZ, 3 }, "qpoints", (char*[]) { "nq", "pol" },
              NULL);
    if ((nc_err = nc_put_var(ncid, varid, phonon->qpts_BZ)))
    {
        ERR(nc_err);
    }

    // qpoints in iBZ
    def_ncVar(ncid, &varid, 2, ELPH_NC4_IO_FLOAT,
              (ND_int[]) { phonon->nq_iBZ, 3 }, "qpoints_iBZ", (char*[]) { "nq_iBZ", "pol" },
              NULL);
    if ((nc_err = nc_put_var(ncid, varid, phonon->qpts_iBZ)))
    {
        ERR(nc_err);
    }

    // write qmap for qpoints (for each qpt in BZ, it gives corresponding iBZ and symm used)
    def_ncVar(ncid, &varid, 2, NC_INT,
              (ND_int[]) { phonon->nq_BZ, 2 }, "qmap", (char*[]) { "nq", "dim_two" },
              NULL);
    if ((nc_err = nc_put_var(ncid, varid, phonon->qmap)))
    {
        ERR(nc_err);
    }

    // write kmap for kpoints (for each kpt in BZ, it gives corresponding iBZ and symm used)
    // This is internally available in yambo, it can be used to cross check if rotation is
    // done in same way

    def_ncVar(ncid, &varid, 2, NC_INT,
              (ND_int[]) { lattice->nkpts_BZ, 2 }, "kmap", (char*[]) { "nk", "dim_two" },
              NULL);
    if ((nc_err = nc_put_var(ncid, varid, lattice->kmap)))
    {
        ERR(nc_err);
    }

    // write start and end band indices
    int bands_tmp[2] = { lattice->start_band, lattice->end_band };
    def_ncVar(ncid, &varid, 1, NC_INT,
              (ND_int[]) { 2 }, "bands", (char*[]) { "two_scalars" },
              NULL);
    if ((nc_err = nc_put_var(ncid, varid, bands_tmp)))
    {
        ERR(nc_err);
    }

    int nph_sym = phonon->nph_sym;
    // write number of phonon symmetries
    def_ncVar(ncid, &varid, 1, NC_INT,
              (ND_int[]) { 1 }, "number_of_phonon_symmetries", (char*[]) { "scalar" },
              NULL);
    if ((nc_err = nc_put_var(ncid, varid, &nph_sym)))
    {
        ERR(nc_err);
    }

    int time_rev_present = 0;
    ELPH_float* symm_mats = malloc(sizeof(ELPH_float) * 9 * nph_sym);
    CHECK_ALLOC(symm_mats);

    ELPH_float* tau_vecs = malloc(sizeof(ELPH_float) * 3 * nph_sym);
    CHECK_ALLOC(tau_vecs);

    for (int isym = 0; isym < nph_sym; ++isym)
    {
        memcpy(symm_mats + isym * 9, phonon->ph_syms[isym].Rmat, sizeof(ELPH_float) * 9);
        memcpy(tau_vecs + isym * 3, phonon->ph_syms[isym].tau, sizeof(ELPH_float) * 3);
        if (phonon->ph_syms[isym].time_rev)
        {
            time_rev_present = 1;
        }
    }

    // write information about the time reversal symmetry of the phonon
    def_ncVar(ncid, &varid, 1, NC_INT,
              (ND_int[]) { 1 }, "time_reversal_phonon", (char*[]) { "scalar" },
              NULL);
    if ((nc_err = nc_put_var(ncid, varid, &time_rev_present)))
    {
        ERR(nc_err);
    }

    // write info about the symmetry_matrices in cart units
    def_ncVar(ncid, &varid, 3, ELPH_NC4_IO_FLOAT,
              (ND_int[]) { nph_sym, 3, 3 }, "symmetry_matrices", (char*[]) { "nsym_ph", "pol", "pol" },
              NULL);
    if ((nc_err = nc_put_var(ncid, varid, symm_mats)))
    {
        ERR(nc_err);
    }

    // write info about the fractional_translation in cart units
    def_ncVar(ncid, &varid, 2, ELPH_NC4_IO_FLOAT,
              (ND_int[]) { nph_sym, 3 }, "fractional_translation", (char*[]) { "nsym_ph", "pol" },
              NULL);

    if ((nc_err = nc_put_var(ncid, varid, tau_vecs)))
    {
        ERR(nc_err);
    }

    // write what type of kernel is used in the calculation
    size_t str_size_tmp = strlen(kernel_str) + 1;
    // for strings we also write the null terminator to the netcdf variable
    def_ncVar(ncid, &varid, 1, NC_CHAR,
              (ND_int[]) { str_size_tmp }, "kernel", (char*[]) { "kernel_str_size" },
              NULL);
    if ((nc_err = nc_put_var(ncid, varid, kernel_str)))
    {
        ERR(nc_err);
    }

    // write the information about the convention used in the code
    str_size_tmp = strlen(convention_str) + 1;
    def_ncVar(ncid, &varid, 1, NC_CHAR,
              (ND_int[]) { str_size_tmp }, "convention", (char*[]) { "convention_str_size" },
              NULL);

    if ((nc_err = nc_put_var(ncid, varid, convention_str)))
    {
        ERR(nc_err);
    }

    free(symm_mats);
    free(tau_vecs);
}
