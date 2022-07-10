#ifndef CORE_H
#define CORE_H

#include <tfhe/tfhe.h>

typedef LweSample LS;
typedef TFheGateBootstrappingCloudKeySet CK;

struct vcf {
    int pos;
    int reflen;
    int altlen;
    int alt[96];
};

struct cvcf {
    LS *pos;
    LS *reflen;
    LS *altlen;
    LS *alt;
};

#endif