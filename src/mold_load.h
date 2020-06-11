#ifndef _MOLD_PARSE_H_
#define _MOLD_PARSE_H_

#include "mold.h"

#ifdef __cplusplus
extern "C" {
#endif

mold_molecule* parse_pdb(const char* str, uint32_t len);
mold_molecule* parse_gro(const char* str, uint32_t len);

#ifdef __cplusplus
}
#endif

#endif