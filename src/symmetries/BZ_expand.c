#include "symmetries.h"

static ELPH_float Function(norm_vec, Nd_cmplxS)(ELPH_float* vec);
static bool Function(check_zero, Nd_cmplxS)(ELPH_float* vec1, ELPH_float* vec2,
                                            ELPH_float r_tol);
static bool Function(Unique_element_finder,
                     Nd_cmplxS)(ELPH_float* array_in, ELPH_float* check_element,
                                ND_int size_array_in);
static void Function(rec_vecs, Nd_cmplxS)(const ELPH_float* lat_vec,
                                          ELPH_float* rec_vec);
static void Function(crystal2cart, Nd_cmplxS)(const ELPH_float* rec_vec,
                                              const ELPH_float* kcrys,
                                              ELPH_float* kcart);
static void Function(cart2crystal, Nd_cmplxS)(const ELPH_float* lat_vec,
                                              const ELPH_float* kcart,
                                              ELPH_float* kcrys);

void bz_expand(ND_array(Nd_floatS) * ibz_kpts, ND_array(Nd_floatS) * sym_mats,
               ND_array(Nd_floatS) * lat_vec, ND_array(Nd_floatS) * kpoints,
               nd_arr_i* kmap)
{
    /*
    **
    Note: It expects the Number of kpoints in full BZ is already known
    INPUTS:
    ibz_kpts --> All k-points in irreducible BZ in cartisian coordinates (3,
    nibz) sym_mats --> All symmetric matrices in cartisian coordinates
    (nsym,3,3) (note these are transposed) a_lat    --> lattice vectors (3,3).
    a[i] = a_lat[:,i]

    OUTPUTS:
    kpoints --> in crystal coorddinates (expanded in full BZ) (nBZ, 3)
    kmap    --> (nbz,2) MAP array. (nbz,0) is related kpoint in iBZ and
                (nbz,1) symmetry mat number that used to transform this.

    k' = Sk (in cart coord)
    k'_red = b^{-1}@k' = (a_lat^T@S)@k where a_lat[:,i] is ith lattice vector
    */

    /* 1. Computing k.T@S.T@alat --> (ibz,nsym,3) */

    ELPH_float blat[9];

    Function(rec_vecs, Nd_cmplxS)(lat_vec->data,
                                  blat);  // get reciprocal vectors

    ND_int Nbz = kpoints->dims[0];
    ND_int Nsym = sym_mats->dims[0];
    ND_int Nibz = ibz_kpts->dims[1];

    ND_array(Nd_floatS) kdotS, krotated;

    ND_function(init, Nd_floatS)(&kdotS, 2, nd_idx{Nibz, 3});
    ND_function(init, Nd_floatS)(&krotated, 3, nd_idx{Nsym, Nibz, 3});

    ND_function(calloc, Nd_floatS)(&kdotS);
    ND_function(calloc, Nd_floatS)(&krotated);

    for (ND_int i = 0; i < Nsym; ++i)
    {
        /* k.T@S*/
        ND_function(matmul, Nd_floatS)('T', 'N', ibz_kpts, sym_mats, &kdotS,
                                       1.0, 0.0, NULL, nd_idx{i}, NULL);
        /* (k.T@S) @ lat*/
        ND_function(matmul, Nd_floatS)('N', 'N', &kdotS, lat_vec, &krotated,
                                       1.0, 0.0, NULL, NULL, nd_idx{i});
    }

    /* 2. Get all the unique k-points */

    ND_int kcounter = 0;

    for (ND_int i = 0; i < Nibz; ++i)
    {
        for (ND_int j = 0; j < Nsym; ++j)
        {
            ELPH_float* k_temp =
                ND_function(ele, Nd_floatS)(&krotated, nd_idx{j, i, 0});

            ELPH_float k_cart[3];  // get cartisian version of k_temp;

            Function(crystal2cart, Nd_cmplxS)(blat, k_temp, k_cart);

            ELPH_float norm = Function(norm_vec, Nd_cmplxS)(k_cart);

            /* Minimize the norm (the bring it back to 1st BZ)*/
            for (int ni = -2; ni < 3; ++ni)
            {
                for (int nj = -2; nj < 3; ++nj)
                {
                    for (int nk = -2; nk < 3; ++nk)
                    {
                        ELPH_float k_buffer[3], k_bufferCC[3];

                        k_buffer[0] = k_temp[0] - ni;
                        k_buffer[1] = k_temp[1] - nj;
                        k_buffer[2] = k_temp[2] - nk;

                        Function(crystal2cart, Nd_cmplxS)(blat, k_buffer,
                                                          k_bufferCC);

                        ELPH_float norm_buffer =
                            Function(norm_vec, Nd_cmplxS)(k_bufferCC);

                        if (norm_buffer < norm - 1E-5)
                        {
                            norm = norm_buffer;
                            memcpy(k_temp, k_buffer, sizeof(ELPH_float) * 3);
                        }
                    }
                }
            }

            if (Function(Unique_element_finder, Nd_cmplxS)(kpoints->data,
                                                           k_temp, kcounter))
            {
                memcpy(kpoints->data + kcounter * 3, k_temp,
                       3 * sizeof(ELPH_float));

                *nd_ele_i(kmap, nd_idx{kcounter, 0}) = i;

                *nd_ele_i(kmap, nd_idx{kcounter, 1}) = j;

                ++kcounter;
            }
        }
    }

    if (kcounter != Nbz)
    {
        printf(" K-point expansion failed \n");
        exit(EXIT_FAILURE);
    }

    ND_function(free, Nd_floatS)(&kdotS);
    ND_function(free, Nd_floatS)(&krotated);
    ND_function(uninit, Nd_floatS)(&kdotS);
    ND_function(uninit, Nd_floatS)(&krotated);
}

/* This Function finds the K+Q indices in kpoint grid */
void get_KplusQ_idxs(ND_array(Nd_floatS) * kpoints, int* KplusQidxs,
                     ELPH_float* Q_pt, ND_array(Nd_floatS) * lat_vec,
                     bool Qincrystal)
{
    /* returns index of the k+q point in the kpoints*/
    /* kpoints (nbz,3)  must be in crystal coordinates.
    and KplusQ_pt (3) .
    Q point Can be in cartisian (set Qincrystal = false), if in crystal set it
    to true alat --> lattice vectors (a[:,i])
    */

    if (!Qincrystal)
    {
        ELPH_float Q_buffer[3];

        Function(cart2crystal, Nd_cmplxS)(lat_vec->data, Q_pt, Q_buffer);

        memcpy(Q_pt, Q_buffer, sizeof(ELPH_float) * 3);
    }

    ELPH_float blat[9];

    Function(rec_vecs, Nd_cmplxS)(lat_vec->data,
                                  blat);  // get reciprocal vectors

    ND_int Nbz = kpoints->dims[0];

    for (ND_int i = 0; i < Nbz; ++i)
    {
        ELPH_float* ktemp;
        ELPH_float KplusQ[3];

        ktemp = ND_function(ele, Nd_floatS)(kpoints, nd_idx{i, 0});

        KplusQ[0] = ktemp[0] + Q_pt[0];
        KplusQ[1] = ktemp[1] + Q_pt[1];
        KplusQ[2] = ktemp[2] + Q_pt[2];

        bool KplusQ_found = false;

        ND_int idx_temp;
        /* Now check the index of K+Q point in the kpoints */
        if (Function(isVECpresent, Nd_cmplxS)(kpoints, KplusQ, &idx_temp))
        {
            KplusQ_found = true;
            KplusQidxs[i] = idx_temp;
        }

        if (!KplusQ_found)
        {
            printf(
                "K+Q point cannot be found in kgrid. The k-grid and q-grid is "
                "not commensurate! \n");
            printf(" The Following is the K+Q point ( %f , %f , %f)! \n",
                   KplusQ[0], KplusQ[1], KplusQ[2]);
            exit(EXIT_FAILURE);
        }
    }
}

bool Function(isVECpresent, Nd_cmplxS)(const ND_array(Nd_floatS) * array,
                                       const ELPH_float* vec, ND_int* idx)
{
    /*
    Checks if kpoint is in given array
    all kpoints should be in crystal coordinates
    */
    bool vec_is_in = false;

    for (ND_int i = 0; i < array->dims[0]; ++i)
    {
        ELPH_float buffer[3];

        ELPH_float* tempj = ND_function(ele, Nd_floatS)(array, nd_idx{i, 0});

        buffer[0] = tempj[0] - vec[0];
        buffer[1] = tempj[1] - vec[1];
        buffer[2] = tempj[2] - vec[2];

        buffer[0] -= rint(buffer[0]);
        buffer[1] -= rint(buffer[1]);
        buffer[2] -= rint(buffer[2]);

        if (Function(norm_vec, Nd_cmplxS)(buffer) < 1E-4)
        {
            *idx = i;
            vec_is_in = true;
            break;
        }
    }

    return vec_is_in;
}

/********************* Helper Functions *********************/

static void Function(rec_vecs, Nd_cmplxS)(const ELPH_float* lat_vec,
                                          ELPH_float* rec_vec)
{
    /*
    Computes reciprocal vectors: from a i.e a[:,i]
    out-put -> b[:,i]
    */

    ELPH_float det =
        +lat_vec[3 * 0 + 0] * (lat_vec[3 * 1 + 1] * lat_vec[3 * 2 + 2] -
                               lat_vec[3 * 2 + 1] * lat_vec[3 * 1 + 2]) -
        lat_vec[3 * 0 + 1] * (lat_vec[3 * 1 + 0] * lat_vec[3 * 2 + 2] -
                              lat_vec[3 * 1 + 2] * lat_vec[3 * 2 + 0]) +
        lat_vec[3 * 0 + 2] * (lat_vec[3 * 1 + 0] * lat_vec[3 * 2 + 1] -
                              lat_vec[3 * 1 + 1] * lat_vec[3 * 2 + 0]);

    if (fabs(det) < 1E-3)
    {
        printf(
            "Inconsistant lattice vectors. Lattice vectors must be linearly "
            "independent ! \n");
        exit(EXIT_FAILURE);
    }

    ELPH_float invdet = 1 / det;

    rec_vec[3 * 0 + 0] = (lat_vec[3 * 1 + 1] * lat_vec[3 * 2 + 2] -
                          lat_vec[3 * 2 + 1] * lat_vec[3 * 1 + 2]) *
                         invdet;
    rec_vec[3 * 1 + 0] = -(lat_vec[3 * 0 + 1] * lat_vec[3 * 2 + 2] -
                           lat_vec[3 * 0 + 2] * lat_vec[3 * 2 + 1]) *
                         invdet;
    rec_vec[3 * 2 + 0] = (lat_vec[3 * 0 + 1] * lat_vec[3 * 1 + 2] -
                          lat_vec[3 * 0 + 2] * lat_vec[3 * 1 + 1]) *
                         invdet;
    rec_vec[3 * 0 + 1] = -(lat_vec[3 * 1 + 0] * lat_vec[3 * 2 + 2] -
                           lat_vec[3 * 1 + 2] * lat_vec[3 * 2 + 0]) *
                         invdet;
    rec_vec[3 * 1 + 1] = (lat_vec[3 * 0 + 0] * lat_vec[3 * 2 + 2] -
                          lat_vec[3 * 0 + 2] * lat_vec[3 * 2 + 0]) *
                         invdet;
    rec_vec[3 * 2 + 1] = -(lat_vec[3 * 0 + 0] * lat_vec[3 * 1 + 2] -
                           lat_vec[3 * 1 + 0] * lat_vec[3 * 0 + 2]) *
                         invdet;
    rec_vec[3 * 0 + 2] = (lat_vec[3 * 1 + 0] * lat_vec[3 * 2 + 1] -
                          lat_vec[3 * 2 + 0] * lat_vec[3 * 1 + 1]) *
                         invdet;
    rec_vec[3 * 1 + 2] = -(lat_vec[3 * 0 + 0] * lat_vec[3 * 2 + 1] -
                           lat_vec[3 * 2 + 0] * lat_vec[3 * 0 + 1]) *
                         invdet;
    rec_vec[3 * 2 + 2] = (lat_vec[3 * 0 + 0] * lat_vec[3 * 1 + 1] -
                          lat_vec[3 * 1 + 0] * lat_vec[3 * 0 + 1]) *
                         invdet;
}

static ELPH_float Function(norm_vec, Nd_cmplxS)(ELPH_float* vec)
{
    // Function to find |vec1|
    return sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
}

static bool Function(check_zero, Nd_cmplxS)(ELPH_float* vec1, ELPH_float* vec2,
                                            ELPH_float r_tol)
{
    // Check if vec1 and vec2 are identical

    ELPH_float vec_temp[3];

    vec_temp[0] = (vec1[0] - vec2[0]) - round(vec1[0] - vec2[0]);
    vec_temp[1] = (vec1[1] - vec2[1]) - round(vec1[1] - vec2[1]);
    vec_temp[2] = (vec1[2] - vec2[2]) - round(vec1[2] - vec2[2]);

    if (Function(norm_vec, Nd_cmplxS)(vec_temp) <= r_tol)
    {
        return true;
    }
    else
    {
        return false;
    }
}

static bool Function(Unique_element_finder,
                     Nd_cmplxS)(ELPH_float* array_in, ELPH_float* check_element,
                                ND_int size_array_in)
{
    // returns true/false if k vector is already in the list
    /*
    Note this is bit different from isVECpresent Function.
    This Functions will only check till size_array_in and works only for
    k-points; isVECpresent will check entire array and gives index too
    */
    bool temp_bool = true;

    ELPH_float tol_temp = 1E-5;

    for (ND_int i = 0; i < size_array_in; ++i)
    {
        // The relative tolarence is set to 10^-4
        if (Function(check_zero, Nd_cmplxS)(array_in + 3 * i, check_element,
                                            tol_temp))
        {
            temp_bool = false;
            break;
        }
    }
    return temp_bool;
}

static void Function(crystal2cart, Nd_cmplxS)(const ELPH_float* rec_vec,
                                              const ELPH_float* kcrys,
                                              ELPH_float* kcart)
{
    /*
    b = rec_vec[:,i]

    kcart = b * kcrys
    */
    for (ND_int i = 0; i < 3; ++i)
    {
        kcart[i] = 0;

        for (ND_int j = 0; j < 3; ++j)
        {
            kcart[i] = kcart[i] + rec_vec[3 * i + j] * kcrys[j];
        }
    }
}

static void Function(cart2crystal, Nd_cmplxS)(const ELPH_float* lat_vec,
                                              const ELPH_float* kcart,
                                              ELPH_float* kcrys)
{
    /*
    a = lat_vec[:,i]

    kcrys = a^T * kcart
    */
    for (ND_int i = 0; i < 3; ++i)
    {
        kcrys[i] = 0;

        for (ND_int j = 0; j < 3; ++j)
        {
            kcrys[i] = kcrys[i] + lat_vec[3 * j + i] * kcart[j];
        }
    }
}
