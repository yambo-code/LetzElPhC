#include "../io.h"

/*
File contains function to read dynamical matrix file (in the old format)
and outputs phonon polarization vectors
*/

// only one cpu calls the routine

#define DYN_READ_BUF_SIZE 10000
#define DYN_FLOAT_BUF_SIZE 200


ND_int read_dyn_qe(const char * dyn_file, struct Lattice * lattice, \
    ELPH_float * restrict qpts, ELPH_cmplx * restrict eigs)
{   
    /*
    // reads all the dynamical matrices in the file
    The function return value is number of dynmats found and read
    the qpoints and eigs are updated according
    */

    ND_int nq_found = 0; // function return value

    char * fgets_err; // error code for fgets

    // First, open the file
    FILE * fp = fopen(dyn_file, "r");
    if (fp == NULL)
    {   
        fprintf(stderr, "Opening file %s failed \n",dyn_file);
        error_msg("Unable to open the dynmat file");
    }

    char * read_buf = malloc(sizeof(char)*DYN_READ_BUF_SIZE);
    ELPH_float * read_fbuf = malloc(sizeof(ELPH_float)*DYN_FLOAT_BUF_SIZE);

    // two comment lines
    fgets(read_buf, DYN_READ_BUF_SIZE, fp);
    fgets(read_buf, DYN_READ_BUF_SIZE, fp);
    fgets(read_buf, DYN_READ_BUF_SIZE, fp); // this line as 9 floats
    if (parser_doubles_from_string(read_buf, read_fbuf) != 9) \
                            error_msg("Error reading line 3 in dyn file");
    ND_int ntype = rint(read_fbuf[0]);
    ND_int natom = rint(read_fbuf[1]);
    ND_int nmodes = natom*3;
    if (natom != lattice->atomic_pos->dims[0]) \
                            error_msg("Wrong number of atoms in dyn file");
    ND_int ibrav = rint(read_fbuf[2]);
    if (ibrav == 0)
    {
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // scratch
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // lat vec 1
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // lat vec 2
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // lat vec 3
    }
    
    ELPH_float * atm_mass = malloc(sizeof(ELPH_float)*(natom+ ntype));
    if(atm_mass == NULL) error_msg("Failed to allocate buffer for atomic masses");

    // read atomic mass for each type
    ELPH_float * atm_mass_type = atm_mass + natom;
    for (ND_int i=0; i<ntype; ++i)
    {
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);
        if (parser_doubles_from_string(read_buf, read_fbuf) != 2) \
                error_msg("Failed to read atomic masses from dyn file");
        atm_mass_type[i] = read_fbuf[1];
        if (fabs(atm_mass_type[i]) < ELPH_EPS) error_msg("Zero masses in dynamical file");
    }
    // now set masses for each atom
    for (ND_int i=0; i<natom; ++i)
    {
        fgets(read_buf, DYN_READ_BUF_SIZE, fp);
        if (parser_doubles_from_string(read_buf, read_fbuf) != 5) \
                error_msg("Failed to read atomic masses from dyn file");
        
        ND_int itype = rint(read_fbuf[1]);
        atm_mass[i] = atm_mass_type[itype-1];
    }

    while(true)
    {
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // empty 
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // content strong
        if (! string_start_with(read_buf, "Dynamical", true)) break; // break the loop if dynmaical not found
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // empty 
        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // qpoint

        ELPH_float * restrict qpt_tmp = qpts + nq_found*3 ;
        
        if (parser_doubles_from_string(read_buf, qpt_tmp) !=3) \
                error_msg("error reading qpoint from dynamat files");

        fgets(read_buf, DYN_READ_BUF_SIZE, fp); // empty 

        ELPH_cmplx * restrict eig = eigs + nq_found*nmodes*nmodes;
        
        for (ND_int ia =0; ia < natom; ++ia)
        {
            for (ND_int ib =0; ib < natom; ++ib)
            {   
                fgets(read_buf, DYN_READ_BUF_SIZE, fp); // read ia and ib
                int itmp, jtmp;
                sscanf(read_buf, "%d  %d", &itmp, &jtmp);
                if (itmp != ia || jtmp != ib) \
                        error_msg("error reading dynamical matrix from dynamat files");
                
                // we must divide by 1/sqrt(Ma*Mb)
                ELPH_float inv_mass_sqtr = sqrt(atm_mass[ia]*atm_mass[ib]); 
                inv_mass_sqtr = 1.0/inv_mass_sqtr;

                for (int ix = 0; ix < 3; ++ix)
                {
                    fgets(read_buf, DYN_READ_BUF_SIZE, fp); // read dynmat
                    // get the values of dynamical matrix
                    if (parser_doubles_from_string(read_buf, read_fbuf) != 6) \
                        error_msg("error reading dynamical matrix from dynamat files");
                    for (int i =0 ; i<6; ++i) read_fbuf[i] *= inv_mass_sqtr;
                    for (int iy = 0; iy < 3; ++iy) // // (ia,ix,ib,iy)
                    {   
                        eig[iy + ib*3 + ix*nmodes + ia*3*nmodes] = read_fbuf[2*iy] + I*read_fbuf[2*iy+1];
                    }
                }
            }
        }
        
        // diagonalize the dynamical matrix and divide with sqrt of masses
        

        // update the counter
        ++nq_found ; 
    }

    

    free(atm_mass);
    free(read_fbuf);
    free(read_buf);
    // close the file
    fclose(fp);

    if (nq_found == 0) error_msg("No dynamical matrices found in the dyn file");

    return  nq_found;
}








