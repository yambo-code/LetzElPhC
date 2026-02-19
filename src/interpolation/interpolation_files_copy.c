#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

#include "common/ELPH_POSIX_func.h"
#include "common/ELPH_copy.h"
#include "common/cwalk/cwalk.h"
#include "common/error.h"
#include "interpolation_utilities.h"
#include "io/ezxml/ezxml.h"

#define BUF_SIZE 2048

void copy_ph_save_to_ph_interpolation_qe(const char* ph_save_dir,
                                         const char* ph_interp_dir)
{
    char src_path[BUF_SIZE];
    char dest_path[BUF_SIZE];

    // Create destination directory
    if (ELPH_mkdir(ph_interp_dir) != 0)
    {
        if (errno != EEXIST)
        {
            error_msg("Failed to create ph_interpolation directory.");
        }
    }

    // 1. Copy data-file-schema.xml
    cwk_path_join_multiple(
        (const char*[]){ph_save_dir, "data-file-schema.xml", NULL}, src_path,
        BUF_SIZE);
    cwk_path_join_multiple(
        (const char*[]){ph_interp_dir, "data-file-schema.xml", NULL}, dest_path,
        BUF_SIZE);

    if (copy_files(src_path, dest_path) != 0)
    {
        fprintf(stderr, "Warning: Error copying file %s to %s.\n", src_path,
                dest_path);
    }

    // 2. Parse data-file-schema.xml to find pseudopotentials
    FILE* fp_xml = fopen(src_path, "r");
    if (fp_xml)
    {
        ezxml_t qexml = ezxml_parse_fp(fp_xml);
        if (qexml)
        {
            ezxml_t atom_specs =
                ezxml_get(qexml, "input", 0, "atomic_species", -1);
            if (atom_specs)
            {
                const char* ntyp_str = ezxml_attr(atom_specs, "ntyp");
                if (ntyp_str)
                {
                    int ntype = atoi(ntyp_str);
                    for (int itype = 0; itype < ntype; ++itype)
                    {
                        ezxml_t pseudo_file_xml = ezxml_get(
                            atom_specs, "species", itype, "pseudo_file", -1);
                        if (pseudo_file_xml)
                        {
                            const char* pseudo_filename = pseudo_file_xml->txt;
                            // Copy pseudo file
                            cwk_path_join_multiple(
                                (const char*[]){ph_save_dir, pseudo_filename,
                                                NULL},
                                src_path, BUF_SIZE);
                            cwk_path_join_multiple(
                                (const char*[]){ph_interp_dir, pseudo_filename,
                                                NULL},
                                dest_path, BUF_SIZE);

                            if (copy_files(src_path, dest_path) != 0)
                            {
                                fprintf(
                                    stderr,
                                    "Warning: Error copying pseudopotential "
                                    "file "
                                    "%s.\n",
                                    src_path);
                            }
                        }
                    }
                }
            }
            ezxml_free(qexml);
        }
        else
        {
            fprintf(stderr, "Warning: Error parsing data-file-schema.xml.\n");
        }
        fclose(fp_xml);
    }
    else
    {
        fprintf(stderr,
                "Warning: Could not open %s for parsing pseudopotentials.\n",
                src_path);
    }

    // 3. Copy tensors.xml
    cwk_path_join_multiple((const char*[]){ph_save_dir, "tensors.xml", NULL},
                           src_path, BUF_SIZE);
    cwk_path_join_multiple((const char*[]){ph_interp_dir, "tensors.xml", NULL},
                           dest_path, BUF_SIZE);
    if (copy_files(src_path, dest_path) != 0)
    {
        // tensors.xml might not exist, ignoring error
        fprintf(stdout, "Tensors.xml not found. Skipping...\n");
    }

    // 4. Copy quadrupole.fmt (if exists)
    cwk_path_join_multiple((const char*[]){ph_save_dir, "quadrupole.fmt", NULL},
                           src_path, BUF_SIZE);
    cwk_path_join_multiple(
        (const char*[]){ph_interp_dir, "quadrupole.fmt", NULL}, dest_path,
        BUF_SIZE);
    copy_files(src_path, dest_path);  // Ignoring return value as it is optional

    return;
}
