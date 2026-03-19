#pragma once
#include "dtypes.h"
#include "elphC.h"

void free_lattice_type(struct Lattice* lattice);
void free_phonon_type(struct Phonon* phonon);
void free_Vloc_table_type(struct Vloc_table* vloc_tab);
void free_local_pseudo_type(struct local_pseudo* lpseudo);
void free_Pseudo_type(struct Pseudo* pseudo);
void free_wfc_type(struct WFC* wfc);
