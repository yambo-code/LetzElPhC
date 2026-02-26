// this file contains posix functions. In the entire codebase,
// only here we define function which are not in strict C99 standard
// Its header defintions should go in ELPH_POSIX_func.h
//
#include "ELPH_POSIX_func.h"

#include <sys/stat.h>
#include <sys/types.h>

#include "error.h"
#if defined(_WIN32)
#include <direct.h>
#define mkdir(path, mode) _mkdir(path)
#endif

int ELPH_mkdir(const char* path) { return mkdir(path, 0777); }

int check_dir_exists(const char* dir_path)
{
    // returns 0 if directory exists.
    struct stat dir_info;

    if (stat(dir_path, &dir_info))
    {
        return ERR_DIR_DOES_NOT_EXIST;
    }
    else if (dir_info.st_mode & S_IFDIR)
    {
        return 0;  // dir exists
    }
    else
    {
        return ERR_NOT_A_DIRECTORY;
    }
}
