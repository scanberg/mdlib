/**************************************************************************
**
**  svd3
**
** Quick singular value decomposition as described by:
** A. McAdams, A. Selle, R. Tamstorf, J. Teran and E. Sifakis,
** "Computing the Singular Value Decomposition of 3x3 matrices
** with minimal branching and elementary floating point operations",
**  University of Wisconsin - Madison technical report TR1690, May 2011
**
**	OPTIMIZED CPU VERSION
** 	Implementation by: Eric Jang
**
**  13 Apr 2014
**
**
**  Modified to use some SSE2 intrinsics
**  Minimal changes to support C through pointers and not references
**  Changed interface to array types
**  Robin Skånberg
**  4 Aug 2021
**************************************************************************/


#define _gamma 5.828427124746190097  // FOUR_GAMMA_SQUARED = sqrt(8)+3;
#define _cstar 0.923879532511286756  // cos(pi/8)
#define _sstar 0.382683432365089771  // sin(p/8)
#define EPSILON_FLT 1e-6

#include "svd3.h"

#include <stdbool.h>
#include <xmmintrin.h>

#include <math.h>

static inline float rsqrt(float x) {
    __m128 res = _mm_rsqrt_ss(_mm_set_ps1(x));
    _mm_store1_ps(&x, res);
    return x;
}

inline float accurateSqrt(float x) {
    return sqrtf(x);
}

static inline void condSwap(bool c, float *X, float *Y) {
    // used in step 2
    float Z = *X;
    *X = c ? *Y : *X;
    *Y = c ?  Z : *Y;
}

static inline void condNegSwap(bool c, float *X, float *Y) {
    // used in step 2 and 3
    float Z = -*X;
    *X = c ? *Y : *X;
    *Y = c ?  Z : *Y;
}

// matrix multiplication C = A * B
static inline void multAB(const float A[3][3],
                   //
                   const float B[3][3],
                   //
                   float C[3][3]) {

    C[0][0] = A[0][0] * B[0][0] + A[0][1] * B[1][0] + A[0][2] * B[2][0];
    C[0][1] = A[0][0] * B[0][1] + A[0][1] * B[1][1] + A[0][2] * B[2][1];
    C[0][2] = A[0][0] * B[0][2] + A[0][1] * B[1][2] + A[0][2] * B[2][2];
    C[1][0] = A[1][0] * B[0][0] + A[1][1] * B[1][0] + A[1][2] * B[2][0];
    C[1][1] = A[1][0] * B[0][1] + A[1][1] * B[1][1] + A[1][2] * B[2][1];
    C[1][2] = A[1][0] * B[0][2] + A[1][1] * B[1][2] + A[1][2] * B[2][2];
    C[2][0] = A[2][0] * B[0][0] + A[2][1] * B[1][0] + A[2][2] * B[2][0];
    C[2][1] = A[2][0] * B[0][1] + A[2][1] * B[1][1] + A[2][2] * B[2][1];
    C[2][2] = A[2][0] * B[0][2] + A[2][1] * B[1][2] + A[2][2] * B[2][2];
}

// matrix multiplication C = Transpose[A] * B
static inline void multAtB(const float A[3][3],
                    //
                    const float B[3][3],
                    //
                    float C[3][3]) {
    C[0][0] = A[0][0] * B[0][0] + A[1][0] * B[1][0] + A[2][0] * B[2][0];
    C[0][1] = A[0][0] * B[0][1] + A[1][0] * B[1][1] + A[2][0] * B[2][1];
    C[0][2] = A[0][0] * B[0][2] + A[1][0] * B[1][2] + A[2][0] * B[2][2];
    C[1][0] = A[0][1] * B[0][0] + A[1][1] * B[1][0] + A[2][1] * B[2][0];
    C[1][1] = A[0][1] * B[0][1] + A[1][1] * B[1][1] + A[2][1] * B[2][1];
    C[1][2] = A[0][1] * B[0][2] + A[1][1] * B[1][2] + A[2][1] * B[2][2];
    C[2][0] = A[0][2] * B[0][0] + A[1][2] * B[1][0] + A[2][2] * B[2][0];
    C[2][1] = A[0][2] * B[0][1] + A[1][2] * B[1][1] + A[2][2] * B[2][1];
    C[2][2] = A[0][2] * B[0][2] + A[1][2] * B[1][2] + A[2][2] * B[2][2];
}

static inline void quatToMat3(const float q[4], float M[3][3]) {
    float x = q[0];
    float y = q[1];
    float z = q[2];
    float w = q[3];

    float qxx = x * x;
    float qyy = y * y;
    float qzz = z * z;
    float qxz = x * z;
    float qxy = x * y;
    float qyz = y * z;
    float qwx = w * x;
    float qwy = w * y;
    float qwz = w * z;

    M[0][0] = 1 - 2 * (qyy + qzz);
    M[0][1] = 2 * (qxy - qwz);
    M[0][2] = 2 * (qxz + qwy);
    M[1][0] = 2 * (qxy + qwz);
    M[1][1] = 1 - 2 * (qxx + qzz);
    M[1][2] = 2 * (qyz - qwx);
    M[2][0] = 2 * (qxz - qwy);
    M[2][1] = 2 * (qyz + qwx);
    M[2][2] = 1 - 2 * (qxx + qyy);
}

static inline void approximateGivensQuaternion(float a11, float a12, float a22, float *ch, float *sh) {
    /*
     * Given givens angle computed by approximateGivensAngles,
     * compute the corresponding rotation quaternion.
     */
    *ch = 2 * (a11 - a22);
    *sh = a12;
    bool b = _gamma * *sh * *sh < *ch * *ch;
    float w = rsqrt(*ch * *ch + *sh * *sh);
    *ch = b ? w * *ch : (float)_cstar;
    *sh = b ? w * *sh : (float)_sstar;
}

static inline void jacobiConjugation(const int x, const int y, const int z, float S[3][3], float q[4]) {
    float ch, sh;
    approximateGivensQuaternion(S[0][0], S[1][0], S[1][1], &ch, &sh);

    float scale = ch * ch + sh * sh;
    float a = (ch * ch - sh * sh) / scale;
    float b = (2 * sh * ch) / scale;

    // make temp copy of S
    float _S[3][3];
    _S[0][0] = S[0][0];
    _S[1][0] = S[1][0];
    _S[1][1] = S[1][1];
    _S[2][0] = S[2][0];
    _S[2][1] = S[2][1];
    _S[2][2] = S[2][2];

    // perform conjugation S = Q'*S*Q
    // Q already implicitly solved from a, b
    S[0][0] =  a * ( a * _S[0][0] + b * _S[1][0]) + b * ( a * _S[1][0] + b * _S[1][1]);
    S[1][0] =  a * (-b * _S[0][0] + a * _S[1][0]) + b * (-b * _S[1][0] + a * _S[1][1]);
    S[1][1] = -b * (-b * _S[0][0] + a * _S[1][0]) + a * (-b * _S[1][0] + a * _S[1][1]);
    S[2][0] =  a * _S[2][0] + b * _S[2][1];
    S[2][1] = -b * _S[2][0] + a * _S[2][1];
    S[2][2] = _S[2][2];

    // update cumulative rotation qV
    float tmp[3];
    tmp[0] = q[0] * sh;
    tmp[1] = q[1] * sh;
    tmp[2] = q[2] * sh;
    sh *= q[3];

    q[0] *= ch;
    q[1] *= ch;
    q[2] *= ch;
    q[3] *= ch;

    // (x,y,z) corresponds to ((0,1,2),(1,2,0),(2,0,1))
    // for (p,q) = ((0,1),(1,2),(0,2))
    q[z] += sh;
    q[3] -= tmp[z];  // w
    q[x] += tmp[y];
    q[y] -= tmp[x];

    // re-arrange matrix for next iteration
    _S[0][0] = S[1][1];
    _S[1][0] = S[2][1];
    _S[1][1] = S[2][2];
    _S[2][0] = S[1][0];
    _S[2][1] = S[2][0];
    _S[2][2] = S[0][0];
    S[0][0] = _S[0][0];
    S[1][0] = _S[1][0];
    S[1][1] = _S[1][1];
    S[2][0] = _S[2][0];
    S[2][1] = _S[2][1];
    S[2][2] = _S[2][2];
}

static inline float dist2(float x, float y, float z) { return x * x + y * y + z * z; }

// finds transformation that diagonalizes a symmetric matrix
static inline void jacobiEigenanlysis( 
        // symmetric matrix
        float S[3][3],
        // quaternion representation of V
        float q[4]) {
    q[0] = 0;
    q[1] = 0;
    q[2] = 0;
    q[3] = 1;

    for (int i = 0; i < 4; i++) {
        // we wish to eliminate the maximum off-diagonal element
        // on every iteration, but cycling over all 3 possible rotations
        // in fixed order (p,q) = (1,2) , (2,3), (1,3) still retains
        //  asymptotic convergence
        jacobiConjugation(0, 1, 2, S, q);  // p,q = 0,1
        jacobiConjugation(1, 2, 0, S, q);  // p,q = 1,2
        jacobiConjugation(2, 0, 1, S, q);  // p,q = 0,2
    }
}

static inline void sortSingularValues(  // matrix that we want to decompose
        float B[3][3],
        // sort V simultaneously
        float V[3][3]) {
    float rho1 = dist2(B[0][0], B[1][0], B[2][0]);
    float rho2 = dist2(B[0][1], B[1][1], B[2][1]);
    float rho3 = dist2(B[0][2], B[1][2], B[2][2]);
    bool c;
    c = rho1 < rho2;
    condNegSwap(c, &B[0][0], &B[0][1]);
    condNegSwap(c, &V[0][0], &V[0][1]);
    condNegSwap(c, &B[1][0], &B[1][1]);
    condNegSwap(c, &V[1][0], &V[1][1]);
    condNegSwap(c, &B[2][0], &B[2][1]);
    condNegSwap(c, &V[2][0], &V[2][1]);
    condSwap(c, &rho1, &rho2);
    c = rho1 < rho3;
    condNegSwap(c, &B[0][0], &B[0][2]);
    condNegSwap(c, &V[0][0], &V[0][2]);
    condNegSwap(c, &B[1][0], &B[1][2]);
    condNegSwap(c, &V[1][0], &V[1][2]);
    condNegSwap(c, &B[2][0], &B[2][2]);
    condNegSwap(c, &V[2][0], &V[2][2]);
    condSwap(c, &rho1, &rho3);
    c = rho2 < rho3;
    condNegSwap(c, &B[0][1], &B[0][2]);
    condNegSwap(c, &V[0][1], &V[0][2]);
    condNegSwap(c, &B[1][1], &B[1][2]);
    condNegSwap(c, &V[1][1], &V[1][2]);
    condNegSwap(c, &B[2][1], &B[2][2]);
    condNegSwap(c, &V[2][1], &V[2][2]);
}

static inline void QRGivensQuaternion(float a1, float a2, float *ch, float *sh) {
    // a1 = pivot point on diagonal
    // a2 = lower triangular entry we want to annihilate
    float epsilon = (float)EPSILON_FLT;
    float rho = accurateSqrt(a1 * a1 + a2 * a2);

    *sh = rho > epsilon ? a2 : 0;
    *ch = fabsf(a1) + fmaxf(rho, epsilon);
    bool b = a1 < 0;
    condSwap(b, sh, ch);
    float w = rsqrt(*ch * *ch + *sh * *sh);
    *ch *= w;
    *sh *= w;
}

static inline void QRDecomposition( 
        // matrix that we want to decompose
        float B[3][3],
        // output Q
        float Q[3][3],
        // output R
        float R[3][3]) {

    float ch1, sh1, ch2, sh2, ch3, sh3;
    float a, b;

    // first givens rotation (ch,0,0,sh)
    QRGivensQuaternion(B[0][0], B[1][0], &ch1, &sh1);
    a = 1 - 2 * sh1 * sh1;
    b = 2 * ch1 * sh1;
    // apply B = Q' * B
    R[0][0] =  a * B[0][0] + b * B[1][0];
    R[0][1] =  a * B[0][1] + b * B[1][1];
    R[0][2] =  a * B[0][2] + b * B[1][2];
    R[1][0] = -b * B[0][0] + a * B[1][0];
    R[1][1] = -b * B[0][1] + a * B[1][1];
    R[1][2] = -b * B[0][2] + a * B[1][2];
    R[2][0] = B[2][0];
    R[2][1] = B[2][1];
    R[2][2] = B[2][2];

    // second givens rotation (ch,0,-sh,0)
    QRGivensQuaternion(R[0][0], R[2][0], &ch2, &sh2);
    a = 1 - 2 * sh2 * sh2;
    b = 2 * ch2 * sh2;
    // apply B = Q' * B;
    B[0][0] = a * R[0][0] + b * R[2][0];
    B[0][1] = a * R[0][1] + b * R[2][1];
    B[0][2] = a * R[0][2] + b * R[2][2];
    B[1][0] = R[1][0];
    B[1][1] = R[1][1];
    B[1][2] = R[1][2];
    B[2][0] = -b * R[0][0] + a * R[2][0];
    B[2][1] = -b * R[0][1] + a * R[2][1];
    B[2][2] = -b * R[0][2] + a * R[2][2];

    // third givens rotation (ch,sh,0,0)
    QRGivensQuaternion(B[1][1], B[2][1], &ch3, &sh3);
    a = 1 - 2 * sh3 * sh3;
    b = 2 * ch3 * sh3;
    // R is now set to desired value
    R[0][0] = B[0][0];
    R[0][1] = B[0][1];
    R[0][2] = B[0][2];
    R[1][0] =  a * B[1][0] + b * B[2][0];
    R[1][1] =  a * B[1][1] + b * B[2][1];
    R[1][2] =  a * B[1][2] + b * B[2][2];
    R[2][0] = -b * B[1][0] + a * B[2][0];
    R[2][1] = -b * B[1][1] + a * B[2][1];
    R[2][2] = -b * B[1][2] + a * B[2][2];

    // construct the cumulative rotation Q=Q1 * Q2 * Q3
    // the number of floating point operations for three quaternion multiplications
    // is more or less comparable to the explicit form of the joined matrix.
    // certainly more memory-efficient!
    float sh12 = sh1 * sh1;
    float sh22 = sh2 * sh2;
    float sh32 = sh3 * sh3;

    Q[0][0] = (-1 + 2 * sh12) * (-1 + 2 * sh22);
    Q[0][1] = 4 * ch2 * ch3 * (-1 + 2 * sh12) * sh2 * sh3 + 2 * ch1 * sh1 * (-1 + 2 * sh32);
    Q[0][2] = 4 * ch1 * ch3 * sh1 * sh3 - 2 * ch2 * (-1 + 2 * sh12) * sh2 * (-1 + 2 * sh32);

    Q[1][0] = 2 * ch1 * sh1 * (1 - 2 * sh22);
    Q[1][1] = -8 * ch1 * ch2 * ch3 * sh1 * sh2 * sh3 + (-1 + 2 * sh12) * (-1 + 2 * sh32);
    Q[1][2] = -2 * ch3 * sh3 + 4 * sh1 * (ch3 * sh1 * sh3 + ch1 * ch2 * sh2 * (-1 + 2 * sh32));

    Q[2][0] = 2 * ch2 * sh2;
    Q[2][1] = 2 * ch3 * (1 - 2 * sh22) * sh3;
    Q[2][2] = (-1 + 2 * sh22) * (-1 + 2 * sh32);
}

void svd(
        // input A
        const float A[3][3],
        // output U
        float U[3][3],
        // output S
        float S[3][3],
        // output V
        float V[3][3]) {

    // normal equations matrix
    float ATA[3][3];

    multAtB(A, A, ATA);

    // symmetric eigenalysis
    float q[4];
    jacobiEigenanlysis(ATA, q);
    quatToMat3(q, V);

    float B[3][3];
    multAB(A, V, B);

    // sort singular values and find V
    sortSingularValues(B, V);

    // QR decomposition
    QRDecomposition(B, U, S);
}

/// polar decomposition can be reconstructed trivially from SVD result
// A = UP
void pd(
        // input A
        const float A[3][3],
        // output U
        float U[3][3],
        // output P
        float P[3][3]) {
    float W[3][3];
    float S[3][3];
    float V[3][3];

    svd(A, W, S, V);

    // P = VSV'
    float T[3][3];
    multAB(V, S, T);

    multAB(T, V, P);

    // U = WV'
    multAB(W, V, U);
}