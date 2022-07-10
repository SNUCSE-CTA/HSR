
#ifndef SEQ_CIRCUIT_H
#define SEQ_CIRCUIT_H

#include <tfhe/tfhe.h>

typedef LweSample LS;

typedef TFheGateBootstrappingCloudKeySet CK;

namespace seq {

void add(LS *c, const LS *a, const LS *b, const int w, const CK *ck);

void compare(LS *c, const LS *a, const LS *b, const int w, const CK *ck);

void enc_const(LS *c, const int a, int w, const CK *ck);

void equals(LS *c, const LS *a, const LS *b, const int w, const CK *ck);

void equals_zero(LS *c, const LS *a, const int w, const CK *ck);

void subtract(LS *c, const LS *a, const LS *b, const int w, const CK *ck);

void subtract_one(LS *c, const LS *a, const int w, const CK *ck);

void multiply_const(LS *c, const LS *a, const int t, const int w, const CK *ck);

void select(LS *c, const LS *s, const LS *a, const LS *b, int w, const CK *ck);
void select_zero(LS *c, const LS *s, const LS *a, int w, const CK *ck);

void copy(LS *a, const LS *b, const int w, const CK *ck);

void fill_zero(LS *a, int w, const CK *ck);

} // namespace seq
#endif