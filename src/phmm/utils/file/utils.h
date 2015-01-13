#ifndef _UTILS_
#define _UTILS_

#include <stdio.h>

#define __MAX_PATH (10000)

#define DATAPATH_ENV_VAR "DATAPATH"
#define LOCAL_DATA_PATH "./"

char* resolve_data_dir();
bool check_file(char* fp);
void validate_file(char* fp);
char* x_fgets(char* buff, int size, FILE* file);
FILE* open_f(char* fp, char* mode);

#endif // _UTILS_

