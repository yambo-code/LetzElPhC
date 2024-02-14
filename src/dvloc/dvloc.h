#pragma once
#include "../wfc/wfc.h"

ELPH_float Vloc_Gspace(ELPH_float* work_arr, const char cutoff,
                       const ELPH_float Gnorm, const ND_int ngrid,
                       const ELPH_float* Vloc_atomic, const ELPH_float* r_grid,
                       const ELPH_float* rab_grid, const ELPH_float Zval,
                       const ELPH_float eta, const ELPH_float volume);

void add_dvscf_qe(ELPH_cmplx* restrict dVscf, const ELPH_cmplx* dVloc,
                  const struct Lattice* lattice);

void elphLocal(const ELPH_float* qpt, struct WFC* wfcs, struct Lattice* lattice,
               int ikq, int ik, int kqsym, int ksym, ELPH_cmplx* dVlocr,
               const struct ELPH_MPI_Comms* Comm, ELPH_cmplx* elph_kq);

void dVlocq(const ELPH_float* qpt, struct Lattice* lattice,
            struct Pseudo* pseudo, const ELPH_cmplx* eigVec, ELPH_cmplx* Vlocr,
            MPI_Comm commK);

void create_vlocg_table(const struct Lattice* lattice, struct Pseudo* pseudo,
                        const struct ELPH_MPI_Comms* Comm);

void free_vlocg_table(struct Vloc_table* vloc_table);
