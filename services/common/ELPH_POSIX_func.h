#pragma once

// Try to avoid any POSIX headers here.

int check_dir_exists(const char* dir_path);

int ELPH_mkdir(const char* path);
