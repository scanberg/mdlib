/* Copyright (c) 2009-2014, Erik Lindahl & David van der Spoel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef XDRFILE_TRR_H
#define XDRFILE_TRR_H

#ifdef __cplusplus
extern "C" {
#endif

#include "xdrfile.h"

/* All functions return exdrOK if succesfull.
 * (error codes defined in xdrfile.h).
 */

/* This function returns the number of atoms in the trr file in *natoms */
int read_trr_natoms(const char* fn, int* natoms);

/* This function returns the number of frames and the number of atoms
 * in the trr file in *natoms and *nframes.
 * It also returns the starting position of each frame as bytes from the beginning of the file
 * in **offsets, which has to be freed manually.
 */
int read_trr_header(const char* fn, int* natoms, unsigned long* nframes, int64_t** offsets);

/* Read one frame of an open trr file. If either of x,v,f,box are
   NULL the arrays will be read from the file but not used.  */
int read_trr(XDRFILE* xd, int natoms, int* step, float* t, float* lambda, matrix box, rvec* x,
             rvec* v, rvec* f, uint8_t* has_prop);

/* Write a frame to trr file */
int write_trr(XDRFILE* xd, int natoms, int step, float t, float lambda, matrix box, rvec* x,
              rvec* v, rvec* f);

/* Minimum TRR header size.
 *  > int(4) magic
 *  > int(4) slen
 *  > string version = uint(4) n + bytes(n)
 *  > 10xint(4) ir_size, e_size, box_size, vir_size, pres_size,
 *               top_size, sym_size, x_size, v_size, f_size
 *  > int(4) natoms
 *  > int(4) step
 *  > int(4) nre
 *  > float(4)/double(8) t
 *  > float(4)/double(8) lamda
 * For an empty version string (n=0) this adds up to 72 bytes.
 * Default version string is "GMX_trn_file" with n=12, so 84 bytes are typical.
 * It can have 8 bytes more if we have double time and lambda.
 */
#define TRR_MIN_HEADER_SIZE 72

/* Flags to signal the existance of box, positions, velocities and forces in a frame */
#define TRR_HAS_BOX 1        /* 0b0001 Box */
#define TRR_HAS_POSITIONS 2  /* 0b0010 Positions */
#define TRR_HAS_VELOCITIES 4 /* 0b0100 Velocities */
#define TRR_HAS_FORCES 8     /* 0b1000 Forces */

#ifdef __cplusplus
}
#endif

#endif
