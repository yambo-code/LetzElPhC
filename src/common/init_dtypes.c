#include "init_dtypes.h"

#include <mpi.h>
#include <string.h>

#include "dtypes.h"
#include "elphC.h"

/*
 * Note : Always set pointers to NULL even after memset as
 * all 0 bits does not need to represent a NULL pointer.
 * Let us leave these optimizations to compilers.
 * Almost all the compilers are good at removing this.
 * but still let's stick to the stardard !
 * */

void init_wfc_type(struct WFC* wfc)
{
    memset(wfc, 0, sizeof(*wfc));

    wfc->wfc = NULL;
    wfc->gvec = NULL;
    wfc->Fk = NULL;
}

void init_lattice_type(struct Lattice* lattice)
{
    memset(lattice, 0, sizeof(*lattice));

    lattice->atomic_pos = NULL;
    lattice->atom_type = NULL;
    lattice->kpt_iredBZ = NULL;
    lattice->kpt_fullBZ = NULL;
    lattice->kpt_fullBZ_crys = NULL;
    lattice->kmap = NULL;
    lattice->syms = NULL;
    return;
}

void init_phonon_type(struct Phonon* phonon)
{
    memset(phonon, 0, sizeof(*phonon));
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

void init_Vloc_table_type(struct Vloc_table* vloc_tab)
{
    memset(vloc_tab, 0, sizeof(*vloc_tab));
    vloc_tab->g_co = NULL;
    vloc_tab->vlocg = NULL;
    vloc_tab->vploc_co = NULL;

    return;
}

void init_local_pseudo_type(struct local_pseudo* lpseudo)
{
    memset(lpseudo, 0, sizeof(*lpseudo));
    lpseudo->Vloc_atomic = NULL;
    lpseudo->r_grid = NULL;
    lpseudo->rab_grid = NULL;

    return;
}

void init_Pseudo_type(struct Pseudo* pseudo)
{
    memset(pseudo, 0, sizeof(*pseudo));
    pseudo->loc_pseudo = NULL;
    pseudo->PP_table = NULL;
    pseudo->Fsign = NULL;
    pseudo->fCoeff = NULL;
    init_Vloc_table_type(pseudo->vloc_table);
    return;
}
