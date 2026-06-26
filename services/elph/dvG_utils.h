#pragma once
/*
 * License-Identifier: GPL
 *
 * Copyright (C) 2025 The Yambo Team
 *
 * Authors (see AUTHORS file for details): RR AM
 *
 * Helpers for gathering the distributed dVscf real-space potential and
 * transforming it to reciprocal space so that the elph_dvG_fill_fn callback
 * can be invoked from commK rank 0.
 */

#include "common/dtypes.h"
#include "elph.h"
#include <mpi.h>

/*
 * Gather z-slabs of dVscf from all commK ranks onto rank 0, perform a
 * forward 3D FFT for each (mode, spin) slice, and return the result.
 *
 * dVscf layout (C row-major, distributed):
 *   (nmodes, nmag, fft_dims[0], fft_dims[1], nfftz_loc)
 *
 * Return value (only meaningful on commK rank 0; NULL on other ranks):
 *   malloc'd buffer of shape (nmodes, nmag, Nx, Ny, Nz) in G-space (complex).
 *   Caller must free().
 */
ELPH_cmplx* gather_dVscf_and_fft(const ELPH_cmplx* dVscf,
                                   const struct Lattice* lattice,
                                   MPI_Comm commK);
