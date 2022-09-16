#ifndef MD_SVD3_H
#define MD_SVD3_H

#ifdef __cplusplus
extern "C" {
#endif

void svd(
    // input A
    const float A[3][3],
    // output U
    float U[3][3],
    // output S
    float S[3][3],
    // output V
    float V[3][3]);

/// polar decomposition can be reconstructed trivially from SVD result
// A = UP
void pd(
    // input A
    const float A[3][3],
    // output U
    float U[3][3],
    // output P
    float P[3][3]);

#ifdef __cplusplus
}
#endif

#endif
