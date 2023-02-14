#include "eval.h"
#include "par_circuit.h"
#include "seq_circuit.h"

using namespace seq;

// Gap penalty function (affine)
void W(LS *r, const LS *l, int g_o, int g_e, int w_int, const CK *ck) {
    LS *tmp = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *go = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *cmp = new_gate_bootstrapping_ciphertext(ck->params);
    enc_const(go, g_o, w_int, ck);

    subtract_one(tmp, l, w_int, ck);
    multiply_const(tmp, tmp, g_e, w_int, ck);
    add(r, tmp, go, w_int, ck);

    equals_zero(cmp, l, w_int, ck);
    bootsNOT(cmp, cmp, ck);
    select_zero(r, cmp, r, w_int, ck);

    delete_gate_bootstrapping_ciphertext(cmp);
    delete_gate_bootstrapping_ciphertext_array(w_int, tmp);
    delete_gate_bootstrapping_ciphertext_array(w_int, go);
}

// HSR_B-VarScore
void var_score(LS *v, LS *lref, cvcf X, cvcf Y, const int s_m, const int s_s,
               const int g_o, const int g_e, int w_int, int w_alt,
               const CK *ck) {

    LS *m = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *g = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *gx = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *gy = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);

    LS *t1 = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *t2 = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *lmax = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *lmin = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *v1 = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *v2 = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);

    LS *c = new_gate_bootstrapping_ciphertext_array(2, ck->params);

    // Eq. (7)
    compare(c, X.reflen, Y.reflen, w_int, ck);
    select(lref, c, X.reflen, Y.reflen, w_int, ck);
    subtract(t1, X.altlen, X.reflen, w_int, ck);
    subtract(t2, Y.altlen, Y.reflen, w_int, ck);

    // Eq. (12) and (13)
    add(gx, lref, t1, w_int, ck);
    add(gy, lref, t2, w_int, ck);

    // Eq. (8) and (9)
    compare(c, t1, t2, w_int, ck);
    select(lmax, c, t1, t2, w_int, ck);
    select(lmin, c, t2, t1, w_int, ck);

    // Eq. (5)
    add(m, lref, lmin, w_int, ck);
    // Eq. (6)
    subtract(g, lmax, lmin, w_int, ck);

    // Checking x.alt == y.alt
    equals(&c[0], X.altlen, Y.altlen, w_int, ck);
    equals(&c[1], X.alt, Y.alt, w_alt, ck);
    bootsAND(c, &c[0], &c[1], ck);

    // Eq. (10)
    multiply_const(lmax, m, s_m, w_int, ck);
    multiply_const(lmin, m, s_s, w_int, ck);
    W(t1, g, g_o, g_e, w_int, ck);
    select(t2, c, lmax, lmin, w_int, ck);
    add(v1, t1, t2, w_int, ck);

    // Eq. (14)
    W(t1, gx, g_o, g_e, w_int, ck);
    W(t2, gy, g_o, g_e, w_int, ck);
    add(v2, t1, t2, w_int, ck);

    // Eq. (15)
    compare(c, v1, v2, w_int, ck);
    select(v, c, v1, v2, w_int, ck);

    // free memory
    delete_gate_bootstrapping_ciphertext_array(w_int, m);
    delete_gate_bootstrapping_ciphertext_array(w_int, g);
    delete_gate_bootstrapping_ciphertext_array(w_int, gx);
    delete_gate_bootstrapping_ciphertext_array(w_int, gy);
    delete_gate_bootstrapping_ciphertext_array(w_int, t1);
    delete_gate_bootstrapping_ciphertext_array(w_int, t2);
    delete_gate_bootstrapping_ciphertext_array(w_int, lmax);
    delete_gate_bootstrapping_ciphertext_array(w_int, lmin);
    delete_gate_bootstrapping_ciphertext_array(w_int, v1);
    delete_gate_bootstrapping_ciphertext_array(w_int, v2);
    delete_gate_bootstrapping_ciphertext_array(2, c);
}

// HSR_B-Score
void score(LS **V, LS **M, LS **S, LS **E, cvcf *X, cvcf *Y, const int size,
           int R_len, const int s_m, const int s_s, const int g_o,
           const int g_e, int w_int, int w_alt, const CK *ck) {

    // first matched region
    enc_const(S[0], 1, w_int, ck);
    subtract_one(E[0], X[0].pos, w_int, ck);
    multiply_const(M[0], E[0], s_m, w_int, ck);

#pragma omp parallel for
    for (int i = 0; i < size; i++) {
        LS *T = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
        LS *lref = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);

        var_score(V[i], lref, X[i], Y[i], s_m, s_s, g_o, g_e, w_int, w_alt, ck);
        add(S[i + 1], lref, X[i].pos, w_int, ck);

        if (i < size - 1) {
            copy(T, X[i + 1].pos, w_int, ck);
        } else {
            enc_const(T, R_len + 1, w_int, ck);
        }

        subtract_one(E[i + 1], T, w_int, ck);
        subtract(lref, T, S[i + 1], w_int, ck);
        multiply_const(M[i + 1], lref, s_m, w_int, ck);

        delete_gate_bootstrapping_ciphertext_array(w_int, T);
        delete_gate_bootstrapping_ciphertext_array(w_int, lref);
    }
}

// HSR_B-Find
void find(LS *m, LS *s, LS *e, LS **V, LS **M, LS **S, LS **E, int size,
          int w_int, const CK *ck) {
    // max subarray
    LS *score = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *start = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *c = new_gate_bootstrapping_ciphertext(ck->params);
    LS *t = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);

    copy(m, M[0], w_int, ck);
    copy(s, S[0], w_int, ck);
    copy(e, E[0], w_int, ck);
    copy(score, m, w_int, ck);
    copy(start, s, w_int, ck);

    for (int i = 0; i < size; i++) {
        par::add(t, score, V[i], w_int, ck);
        bootsNOT(c, &t[w_int - 1], ck);
        par::select_zero(t, c, t, w_int, ck);
        par::add(score, t, M[i + 1], w_int, ck);
        par::select(start, c, start, S[i + 1], w_int, ck);

        par::compare(c, score, m, w_int, ck);
        par::select(m, c, score, m, w_int, ck);
        par::select(s, c, start, s, w_int, ck);
        par::select(e, c, E[i + 1], e, w_int, ck);
    }

    delete_gate_bootstrapping_ciphertext_array(w_int, start);
    delete_gate_bootstrapping_ciphertext_array(w_int, score);
    delete_gate_bootstrapping_ciphertext_array(w_int, t);
    delete_gate_bootstrapping_ciphertext(c);
}
