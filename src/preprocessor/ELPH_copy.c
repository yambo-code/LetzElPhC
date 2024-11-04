// This file contains functions related to coping files and directories.
// Contains os dependent functions.

#include "ELPH_copy.h"

#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
// above two headers are found in both windows and posix systems.

int check_dir_exists(const char *dir_path)
{
    // returns 0 if directory exists.
    struct stat dir_info;

    if (stat(dir_path, &dir_info))
    {
        return 1;
    }
    else if (dir_info.st_mode & S_IFDIR)
    {
        return 0;  // dir exists
    }
    else
    {
        return 1;
    }
}

int copy_files(const char *file_read, const char *file_write)
{
    // returns 0 if copying is sucuessful.
    FILE *frd = fopen(file_read, "rb");
    if (!frd)
    {
        return 1;
    }
    FILE *fwr = fopen(file_write, "wb");
    if (!fwr)
    {
        fclose(frd);
        return 1;
    }

    char read_buf[BUFSIZ];  // BUFSIZ is defined in <stdio.h>
    size_t nbytes;

    int err_code = 0;
    while ((nbytes = fread(read_buf, 1, sizeof(read_buf), frd)) > 0)
    {
        if (fwrite(read_buf, 1, nbytes, fwr) != nbytes)
        {
            err_code = 1;
            break;
        }
    }

    fclose(frd);
    fclose(fwr);

    return err_code;
}
