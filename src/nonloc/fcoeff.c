/*
These routines contain functions related the pseudo potentials.

!! DO NOT typedef structs !!
*/

#include "fcoeff.h"

/** Static declarations*/
static ELPH_cmplx Cmat(int m1, int m2);
static ELPH_float CGCoeff(bool jp1, bool spin_down, int l, int twomj);
static ELPH_cmplx Umat(bool jp1, bool spin_down, int l, int twomj, int m2);

void alloc_and_Compute_f_Coeff(struct Lattice* lattice, struct Pseudo* pseudo)
{
    /*
    This functions allocates memory and computes f coefficients for spinors
    which is given by f^{l,j,\sigma,\sigma'}_{m,m'} =\sum_{m_j}
    P^{l,j,\sigma}_{m_j,m} * conj(P^{l,j,\sigma'}_{m_j,m'})

    where P^{l,j,\sigma}_{m_j,m} = \alpha^{\sigma,l,j}_{m_j} *
    U^{\sigma,l,j}_{m_j,m} Refer Eq.9 of PHYSICAL REVIEW B 71, 115106 s2005d

    Input : PP_table from yambo. This contains l,j values (l*j,ntype,3)
    (l+1,2j,pp_spin), pp_spin = 1 ? Output : (l*j,ntype) array. each element of
    array has a shape (2l+1, 2l+1, nspinor, nspinor)

    Note, we only compute these matrices once.

    also return no of Coeffs allocated.

    Note: donot forget to free the memory
    */

    int nspinor = lattice->nspinor;

    if (nspinor == 1 || !lattice->is_soc_present)
    {
        pseudo->fCoeff = NULL;
        return;
    }

    ND_int nltimesj = pseudo->nltimesj;
    ND_int ntype = pseudo->ntype; // number of atomic types
    ND_int nproj = ntype * nltimesj;

    ELPH_float* PP_table = pseudo->PP_table;
    pseudo->fCoeff = malloc(sizeof(ELPH_cmplx*) * nproj);
    CHECK_ALLOC(pseudo->fCoeff);
    // (ntype,l*j)
    ELPH_cmplx** Coeff = pseudo->fCoeff;

    // initiate to NULL (this is must else U.B when freeing)
    for (ND_int i = 0; i < nproj; ++i)
    {
        Coeff[i] = NULL;
    }

    for (ND_int i = 0; i < nproj; ++i)
    {
        ND_int lidx = i / ntype;
        ND_int itype = i % ntype;
        // (l*j,ntype,3) [ntype*3,3,1]
        int l = rint(PP_table[ntype * 3 * lidx + itype * 3] - 1); // PP_table[lidx,itype,0]
        int j = rint(
            PP_table[ntype * 3 * lidx + itype * 3 + 1]); // Careful, this is 2j

        if (l < 0 || j ==0 )
        {
            continue;
        } // skip fake entries or j ==0 (j =0 if one of the pseudo does not have soc)

        ND_int two_lp1 = 2 * l + 1;
        ND_int coeff_len = two_lp1 * two_lp1 * nspinor * nspinor;
        Coeff[i] = malloc(sizeof(ELPH_cmplx) * coeff_len); //(two_lp1,two_lp1,nspinor,nspinor)
        CHECK_ALLOC(Coeff[i]);
        // zero out the buffer
        for (ND_int inum = 0; inum < coeff_len; ++inum)
        {
            Coeff[i][inum] = 0;
        }

        bool jp1 = true; // This is true when j = l + 1/2 else false

        if (j == 2 * l + 1)
        {
            jp1 = true;
        }
        else if (j == 2 * l - 1)
        {
            jp1 = false;
        }
        else
        {
            error_msg("Inconsistant j value in Pseudo potential.");
        }

        for (ND_int il1 = 0; il1 < 2 * l + 1; ++il1) // m
        {
            for (ND_int il2 = 0; il2 < 2 * l + 1; ++il2) // m'
            {
                for (ND_int is1 = 0; is1 < nspinor; ++is1) // sigma
                {
                    for (ND_int is2 = 0; is2 < nspinor; ++is2) // sigma'
                    {
                        int im1 = il1 - l;
                        int im2 = il2 - l;
                        /* is =0 for spin up and is = 1 for spin down
                        so we can use "is" as boolean for spin down
                        */
                        ELPH_cmplx sum = 0;

                        for (int mj = -j; mj <= j; mj += 2)
                        { /** Careful ! mj has a factor 2*/
                            sum += CGCoeff(jp1, is1, l, mj) * CGCoeff(jp1, is2, l, mj) * Umat(jp1, is1, l, mj, im1) * conj(Umat(jp1, is2, l, mj, im2));
                        }

                        Coeff[i][nspinor * nspinor * (il1 * two_lp1 + il2) + is1 * nspinor + is2] = sum;
                        // Coeff[i][il1,il2,is1,is2] = sum
                    }
                }
            }
        }
    }
    return;
}

void free_f_Coeff(struct Lattice* lattice, struct Pseudo* pseudo)
{
    if (!pseudo->fCoeff)
    {
        return;
    }

    ND_int nltimesj = pseudo->nltimesj;
    ND_int ntype = pseudo->ntype; // number of atomic types
    ND_int nproj = ntype * nltimesj;

    ELPH_cmplx** Coeff = pseudo->fCoeff;

    for (ND_int i = 0; i < nproj; ++i)
    {
        free(Coeff[i]);
    }
    free(Coeff);
    pseudo->fCoeff = NULL;
}

static ELPH_cmplx Umat(bool jp1, bool spin_down, int l, int twomj, int m2)
{
    /* These are the U matrices entering into Eq.6 of PHYSICAL REVIEW B 71,
     * 115106 s2005 */
    if (abs(m2) > l)
    {
        error_msg("m value greater than l");
    }
    if (jp1)
    {
        int m = (twomj - 1) / 2;
        if (!spin_down)
        {
            return Cmat(m, m2);
        }
        else
        {
            return Cmat(m + 1, m2);
        }
    }
    else
    {
        int m = (twomj + 1) / 2;
        if (!spin_down)
        {
            return Cmat(m - 1, m2);
        }
        else
        {
            return Cmat(m, m2);
        }
    }
    return 0.0;
}

/* Function to compute the transformation matrix from real to complex spherical
 * harmonics. */
static ELPH_cmplx Cmat(int m1, int m2)
{
    /*
    This function returns transformation coeffients for transforming
    Real to complex spherical harmonics. These are C^\dagger elements given in
    Eq. 19 of Ref:
    Miguel A. Blanco, M. Flórez, M. Bermejo
    Evaluation of the rotation matrices in the basis of real spherical harmonics
    J. Mol. Struct. THEOCHEM 419, 19-27 (1997)

    It should be noted that definition of these elements are implementation
    specific.
    */
    int abs_m1 = abs(m1);
    int abs_m2 = abs(m2);

    if (abs_m1 != abs_m2)
    {
        return 0.0;
    }
    if (m1 == 0 && m2 == 0)
    {
        return 1.0;
    }
    if (m1 == abs_m1 && m2 == abs_m2)
    {
        return pow(-1.0, abs_m1) / ELPH_SQRT2;
    }
    if (m1 == -abs_m1 && m2 == abs_m2)
    {
        return 1.0 / ELPH_SQRT2;
    }
    if (m1 == abs_m1 && m2 == -abs_m2)
    {
        return I * pow(-1.0, abs_m1) / ELPH_SQRT2;
    }
    if (m1 == -abs_m1 && m2 == -abs_m2)
    {
        return -I / ELPH_SQRT2;
    }
    else
    {
        error_msg("Incorrect m1 and m2 values");
    }
    return 0;
}

static ELPH_float CGCoeff(bool jp1, bool spin_down, int l, int twomj)
{ /*
  Function to compute Clebsch-Gordan coefficients alpha
  jp1 is true for j = l+1/2 else false.
  Eq-5 of Ref : PHYSICAL REVIEW B 71, 115106 s2005d
  */
    /* check if integrand is non negative */
    /* Warning twomj = 2*j, we pass 2j instead of j*/
    ELPH_float mj = twomj / 2.0;
    if (jp1)
    {
        ELPH_float m = mj - 0.5;
        if (!spin_down)
        {
            return sqrt((l + m + 1.0) / (2.0 * l + 1.0));
        }
        else
        {
            return sqrt((l - m) / (2.0 * l + 1.0));
        }
    }
    else
    {
        ELPH_float m = mj + 0.5;
        if (!spin_down)
        {
            return sqrt((l - m + 1.0) / (2.0 * l + 1.0));
        }
        else
        {
            return -sqrt((l + m) / (2.0 * l + 1.0));
        }
    }
}
