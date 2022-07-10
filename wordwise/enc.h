#ifndef ENC_H
#define ENC_H

#include "core.h"
#include <string>

void encrypt(cvcf &X, cvcf &Y, string xfile, string yfile, int size, int w_alt,
             int n, Encryptor enc, PublicKeyPack keypack);

#endif
