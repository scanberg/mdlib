#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Returns the length in bytes of the encoded string which is required to be the minimum buffer size supplied as 'output' into base64_encode
size_t md_base64_encode_size_in_bytes(size_t length);

// Encodes the given input bytes as Base64
// The output is expected to have a capacity as mentioned above.
// Returns the encoded length in bytes
size_t md_base64_encode(char* out_encoded, const void* in_raw, size_t length);

// Gives the length in bytes of the encoded string which is required to be the minimum buffer size supplied as 'output' into base64_decode
size_t md_base64_decode_size_in_bytes(size_t length);

// Decodes a sequence of characters given in Base64
// The output is expected to have a capacity as mentioned above.
// Returns the decoded length in bytes
size_t md_base64_decode(void* out_raw, const char *in_encoded, size_t length);

#ifdef __cplusplus
}
#endif
