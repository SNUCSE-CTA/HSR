#ifndef CORE_H
#define CORE_H
#include "HEAAN.h"
#include "utils/MathUtils.h"

struct vcf {
    int pos;
    int reflen;
    int altlen;
    string alt;
};

struct cvcf {
    Ciphertext pos;
    Ciphertext reflen;
    Ciphertext altlen;
    vector<Ciphertext> alt;
};

#endif
