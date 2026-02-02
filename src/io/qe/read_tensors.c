// reads tensors born charges and dielectric tensor
// from tensors.xml file
//

#include <ctype.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/dtypes.h"
#include "common/error.h"
#include "common/numerical_func.h"
#include "common/string_func.h"
#include "elphC.h"
#include "io/ezxml/ezxml.h"

#define MAX_LINE 1024

void read_ph_tensors_qe(const char* tensor_xml_file, const ND_int natom,
                        struct Phonon* phonon)
{
    phonon->Zborn = NULL;
    phonon->epsilon = NULL;
    phonon->Qpole = NULL;
    // for now no Quadruple as q.e does not support it.

    FILE* fp = fopen(tensor_xml_file, "r");
    if (NULL == fp)
    {
        // file not found or not readable.
        // so we skip without reading them.
        return;
    }

    ezxml_t tensor_xml = ezxml_parse_fp(fp);
    if (NULL == tensor_xml)
    {
        error_msg("parsing tensor.xml file failed\n");
    }

    char* tmp_str;
    // first check if tensors exist !
    //
    bool epsilon_exists = false;
    bool Zeu_exists = false;

    ezxml_t txml;

    txml =
        ezxml_get(tensor_xml, "EF_TENSORS", 0, "DONE_EFFECTIVE_CHARGE_EU", -1);
    if (!txml)
    {
        error_msg(
            "parsing DONE_EFFECTIVE_CHARGE_EU from tensor.xml file failed\n");
    }
    tmp_str = txml->txt;
    if ('t' == tolower((unsigned char)(*tmp_str)))
    {
        Zeu_exists = true;
    }

    txml = ezxml_get(tensor_xml, "EF_TENSORS", 0, "DONE_ELECTRIC_FIELD", -1);
    if (!txml)
    {
        error_msg(
            "parsing DONE_EFFECTIVE_CHARGE_EU from tensor.xml file failed\n");
    }
    tmp_str = txml->txt;
    if ('t' == tolower((unsigned char)(*tmp_str)))
    {
        epsilon_exists = true;
    }

    if (epsilon_exists)
    {
        if (Zeu_exists)
        {
            phonon->Zborn = malloc(9 * natom * sizeof(*phonon->Zborn));
            CHECK_ALLOC(phonon->Zborn);
        }

        phonon->epsilon = malloc(9 * sizeof(*phonon->epsilon));
        CHECK_ALLOC(phonon->epsilon);
        // read dielectric tensor
        txml =
            ezxml_get(tensor_xml, "EF_TENSORS", 0, "DIELECTRIC_CONSTANT", -1);
        if (!txml)
        {
            error_msg(
                "parsing DIELECTRIC_CONSTANT from tensor.xml file failed\n");
        }

        tmp_str = txml->txt;

        if (parse_floats_from_string(tmp_str, phonon->epsilon, 9) != 9)
        {
            error_msg("Parsing epsilon failed");
        }
        // read effective charges
        if (Zeu_exists)
        {
            txml = ezxml_get(tensor_xml, "EF_TENSORS", 0,
                             "EFFECTIVE_CHARGES_EU", -1);
            if (!txml)
            {
                error_msg("parsing Born charges from tensor.xml file failed\n");
            }
            tmp_str = txml->txt;

            if (parse_floats_from_string(tmp_str, phonon->Zborn, 9 * natom) !=
                9 * natom)
            {
                error_msg("Parsing Born charges failed");
            }
        }
        // transpose epsilon and born charges as they are stored in transpose
        // order
        transpose3x3f_inplace(phonon->epsilon);

        if (Zeu_exists)
        {
            for (ND_int ia = 0; ia < natom; ++ia)
            {
                transpose3x3f_inplace(phonon->Zborn + 9 * ia);
            }
        }
    }
    else
    {
        fprintf(stderr,
                "Warning : Tensor file exists but Born charges or dielectric "
                "tensor missing !\n");
    }

    ezxml_free(tensor_xml);
    fclose(fp);
}

void read_quadrupole_fmt(const char* filename, ELPH_float** Qpole_buf,
                         int natom)
{
    // (natom, x,y, atom_dir)
    *Qpole_buf = NULL;
    //
    FILE* fp = fopen(filename, "r");
    if (!fp)
    {
        return;
    }
    //
    char line[MAX_LINE];
    fgets(line, sizeof(line), fp);
    // Read data lines
    ND_int ndata_parsed = 0;
    //
    ELPH_float* buffer = malloc(sizeof(*buffer) * 27 * natom);
    CHECK_ALLOC(buffer);
    *Qpole_buf = buffer;
    //
    while (fgets(line, sizeof(line), fp))
    {
        bool is_empty = true;
        for (ND_int i = 0; line[i] != '\0'; i++)
        {
            if (!isspace((unsigned char)line[i]))
            {
                is_empty = false;
                break;
            }
        }

        // Skip empty lines
        if (is_empty || line[0] == '\n' || line[0] == '\0')
        {
            continue;
        }
        //
        if (ndata_parsed >= 3 * natom)
        {
            free(buffer);
            error_msg("More number of lines given im the file");
        }
        //
        double Qxx, Qyy, Qzz, Qyz, Qxz, Qxy;
        long long atom_idx, dir_idx;

        // Parse the line
        int parsed =
            sscanf(line, "%lld %lld %lf %lf %lf %lf %lf %lf", &atom_idx,
                   &dir_idx, &Qxx, &Qyy, &Qzz, &Qyz, &Qxz, &Qxy);

        --atom_idx;
        --dir_idx;
        if (parsed == 8)
        {
            long long idx;
            // diagonal
            idx = atom_idx * 27 + dir_idx;
            buffer[idx] = Qxx;
            idx = atom_idx * 27 + 1 * 9 + 1 * 3 + dir_idx;
            buffer[idx] = Qyy;
            idx = atom_idx * 27 + 2 * 9 + 2 * 3 + dir_idx;
            buffer[idx] = Qzz;
            // off diagonal
            idx = atom_idx * 27 + 0 * 9 + 1 * 3 + dir_idx;
            buffer[idx] = Qxy;
            idx = atom_idx * 27 + 1 * 9 + 0 * 3 + dir_idx;
            buffer[idx] = Qxy;
            //
            idx = atom_idx * 27 + 0 * 9 + 2 * 3 + dir_idx;
            buffer[idx] = Qxz;
            idx = atom_idx * 27 + 2 * 9 + 0 * 3 + dir_idx;
            buffer[idx] = Qxz;
            // Q[1][2] = Q[2][1] = Qyz
            idx = atom_idx * 27 + 1 * 9 + 2 * 3 + dir_idx;
            buffer[idx] = Qyz;
            idx = atom_idx * 27 + 2 * 9 + 1 * 3 + dir_idx;
            buffer[idx] = Qyz;
        }
        else
        {
            free(buffer);
            error_msg(
                "Unable to parse atom, index Qxx, Qyy, Qzz, Qyz, Qxz, Qxy");
        }
        ++ndata_parsed;
    }
    fclose(fp);
    if (ndata_parsed != 3 * natom)
    {
        free(buffer);
        error_msg("Unable to real full quadrupole tensor.");
    }
}
