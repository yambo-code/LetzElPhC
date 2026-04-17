#pragma once
#include "dtypes.h"
#include "elphC.h"

void init_lattice_type(struct Lattice* lattice);
void init_phonon_type(struct Phonon* phonon);
void init_Vloc_table_type(struct Vloc_table* vloc_tab);
void init_local_pseudo_type(struct local_pseudo* lpseudo);
void init_Pseudo_type(struct Pseudo* pseudo);
void init_wfc_type(struct WFC* wfc);
