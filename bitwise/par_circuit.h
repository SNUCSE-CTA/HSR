
#ifndef PAR_CIRCUIT_H
#define PAR_CIRCUIT_H

#include <tfhe/tfhe.h>

typedef LweSample LS;
typedef TFheGateBootstrappingCloudKeySet CK;

namespace par {

void add(LS *c, const LS *a, const LS *b, const int w, const CK *ck);

void compare(LS *c, const LS *a, const LS *b, const int w, const CK *ck);

void select(LS *c, const LS *s, const LS *a, const LS *b, int w, const CK *ck);

void select_zero(LS *c, const LS *s, const LS *a, int w, const CK *ck);

} // namespace par
#endif