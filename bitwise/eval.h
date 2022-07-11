
#ifndef EVAL_H
#define EVAL_H

#include "core.h"

void score(LS **V, LS **M, LS **S, LS **E, cvcf *X, cvcf *Y, const int size,
           int R_len, const int s_m, const int s_s, const int g_o,
           const int g_e, int w_int, int w_alt, const CK *ck);

void find(LS *m, LS *s, LS *e, LS **V, LS **M, LS **S, LS **E, const int size,
          const int w_int, const CK *ck);

#endif