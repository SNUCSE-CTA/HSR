
#ifndef ENC_H
#define ENC_H

#include "core.h"
#include <string>

void encrypt(cvcf *S, std::string file_name, const int size, const int w,
             const int w_alt, const TFheGateBootstrappingSecretKeySet *key,
             const TFheGateBootstrappingParameterSet *params);

void encrypt_to_file(std::string file_name, const int size, const int w,
                     const int w_alt,
                     const TFheGateBootstrappingSecretKeySet *key,
                     const TFheGateBootstrappingParameterSet *params);
#endif