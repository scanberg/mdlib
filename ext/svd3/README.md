Fast 3x3 SVD
===========

This is an implementation of the method described in <a href="http://pages.cs.wisc.edu/~sifakis/papers/SVD_TR1690.pdf">"Computing the Singular Value Decomposition of 3x3 matrices with minimal branching and elementary floating point operations"</a>. I implemented this as part of <a href="http://wyegelwel.github.io/snow/">a group project</a> for a computer graphics course. 

Execution time per svd call on the CPU is about 2.0 microseconds. Tested on a AMD Phenom(tm) II X4 965 Processor. 

Execution time on the GPU is about 174 microseconds. Tested on a NVIDIA GeForce GTX 460 (profiled using nvvp).

Also included are routines for diagonalization / QR decomposition of 3x3 matrices, which may be useful in their own right. 


##Usage

Just include the header file and you are good to go! 

```C

#include "svd3.h"
float A[3][3];

a[0][0] = -0.558253; a[0][1] = -0.0461681; a[0][2] = -0.505735;
a[1][0] = -0.411397; a[1][1] =  0.0365854; a[1][2] =  0.199707;
a[2][0] =  0.285389; a[2][1] = -0.3137890; a[2][2] =  0.200189;

float U[3][3];
float S[3][3];
float V[3][3];

svd(A, U, S, V);

```

See the included Mathematica notebook for derivations of numerical shortcuts.

## License
MIT License, Eric V. Jang 2014