#pragma once

#ifdef __cplusplus
extern "C" {
#endif

// Returns the length in bytes of the encoded string which is required to be the minimum buffer size supplied as 'output' into base64_encode
int md_base64_encode_size_in_bytes(int length);

// Encodes the given input bytes as Base64
// The output is expected to have a capacity as mentioned above.
// Returns the encoded length in bytes
int md_base64_encode(char* out_encoded_data, const void* in_raw_data, int length);

// Gives the length in bytes of the encoded string which is required to be the minimum buffer size supplied as 'output' into base64_decode
int md_base64_decode_size_in_bytes(int length);

// Decodes a sequence of characters given in Base64
// The output is expected to have a capacity as mentioned above.
// Returns the decoded length in bytes
int md_base64_decode(void* out_raw_data, const char *in_encoded_data, int length);

#ifdef __cplusplus
}
#endif
