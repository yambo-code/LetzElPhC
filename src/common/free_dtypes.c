#include "free_dtypes.h"

#include <mpi.h>
#include <stdlib.h>

#include "dtypes.h"
#include "elphC.h"

void free_wfc_type(struct WFC* Wfc)
{
    free(Wfc->wfc);
    free(Wfc->gvec);
    free(Wfc->Fk);

    Wfc->wfc = NULL;
    Wfc->gvec = NULL;
    Wfc->Fk = NULL;
}

void free_lattice_type(struct Lattice* lattice)
{
    if (!lattice)
    {
        return;
    }

    free(lattice->atomic_pos);
    free(lattice->atom_type);
    free(lattice->kpt_iredBZ);
    free(lattice->kpt_fullBZ);
    free(lattice->kpt_fullBZ_crys);
    free(lattice->kmap);
    free(lattice->syms);

    lattice->atomic_pos = NULL;
    lattice->atom_type = NULL;
    lattice->kpt_iredBZ = NULL;
    lattice->kpt_fullBZ = NULL;
    lattice->kpt_fullBZ_crys = NULL;
    lattice->kmap = NULL;
    lattice->syms = NULL;
    return;
}

void free_phonon_type(struct Phonon* phonon)
{
    if (!phonon)
    {
        return;
    }

    free(phonon->qpts_iBZ);
    free(phonon->qpts_BZ);
    free(phonon->ph_syms);
    free(phonon->qmap);
    free(phonon->nqstar);
    free(phonon->epsilon);
    free(phonon->Zborn);
    free(phonon->Qpole);

    phonon->qpts_iBZ = NULL;
    phonon->qpts_BZ = NULL;
    phonon->ph_syms = NULL;
    phonon->qmap = NULL;
    phonon->nqstar = NULL;
    phonon->epsilon = NULL;
    phonon->Zborn = NULL;
    phonon->Qpole = NULL;

    return;
}

void free_Vloc_table_type(struct Vloc_table* vloc_tab)
{
    if (!vloc_tab)
    {
        return;
    }
    // free the allocate Vloc_table
    free(vloc_tab->g_co);
    free(vloc_tab->vlocg);
    free(vloc_tab->vploc_co);

    vloc_tab->g_co = NULL;
    vloc_tab->vlocg = NULL;
    vloc_tab->vploc_co = NULL;

    return;
}

void free_local_pseudo_type(struct local_pseudo* lpseudo)
{
    if (!lpseudo)
    {
        return;
    }
    free(lpseudo->Vloc_atomic);
    free(lpseudo->r_grid);
    free(lpseudo->rab_grid);
    //
    lpseudo->Vloc_atomic = NULL;
    lpseudo->r_grid = NULL;
    lpseudo->rab_grid = NULL;

    return;
}

void free_Pseudo_type(struct Pseudo* pseudo)
{
    if (!pseudo)
    {
        return;
    }

    // NM: Carefull, pseudo->fCoeff contains
    // data, so the user must free them
    // using free_f_Coeff in fcoeff.c
    // before call this function.
    free(pseudo->fCoeff);
    free_Vloc_table_type(pseudo->vloc_table);
    free(pseudo->Fsign);
    free(pseudo->PP_table);

    if (pseudo->loc_pseudo)
    {
        for (ND_int itype = 0; itype < pseudo->ntype; ++itype)
        {
            free_local_pseudo_type(pseudo->loc_pseudo + itype);
        }
    }
    free(pseudo->loc_pseudo);

    pseudo->loc_pseudo = NULL;
    pseudo->PP_table = NULL;
    pseudo->Fsign = NULL;
    pseudo->fCoeff = NULL;
    return;
}
