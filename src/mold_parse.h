#ifndef _MOLD_PARSE_H_
#define _MOLD_PARSE_H_

#include "mold.h"

#ifdef __cplusplus
extern "C" {
#endif

mold_error mold_parse_gro();
mold_error mold_parse_pdb();

#ifdef __cplusplus
}
#endif

#endif