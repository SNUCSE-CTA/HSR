#ifndef EVAL_H
#define EVAL_H

#include "core.h"

void score(Ciphertext &V, Ciphertext &M, cvcf X, cvcf Y, int n, int R_len,
           int s_m, int s_s, int g_o, int g_e, int w_alt, int d, Ciphertext one,
           HomEvaluator eval, PublicKeyPack keypack);

void find(Ciphertext &score, Ciphertext V, Ciphertext M, int size, int n, int d,
          HomEvaluator eval, PublicKeyPack keypack);

#endif
