#include "ELPH_copy.h"

#include <stdio.h>

#include "common/error.h"

int copy_files(const char* file_read, const char* file_write)
{
    // returns 0 if copying is sucuessful.
    FILE* frd = fopen(file_read, "rb");
    if (!frd)
    {
        return ERR_FILE_OPEN_READ;
    }
    FILE* fwr = fopen(file_write, "wb");
    if (!fwr)
    {
        fclose(frd);
        return ERR_FILE_OPEN_WRITE;
    }

    unsigned char read_buf[BUFSIZ];  // BUFSIZ is defined in <stdio.h>
    size_t nbytes;

    int err_code = 0;
    while ((nbytes = fread(read_buf, 1, sizeof(read_buf), frd)) > 0)
    {
        if (fwrite(read_buf, 1, nbytes, fwr) != nbytes)
        {
            err_code = ERR_FILE_COPY_FAIL;
            break;
        }
    }

    fclose(frd);
    fclose(fwr);

    return err_code;
}
