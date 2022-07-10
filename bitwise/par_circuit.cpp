#include "par_circuit.h"

using namespace par;

// Sklansky matrices for 32-bit adder
const int idx1[5][16] = {
    {1, 3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31},
    {2, 3, 6, 7, 10, 11, 14, 15, 18, 19, 22, 23, 26, 27, 30, 31},
    {4, 5, 6, 7, 12, 13, 14, 15, 20, 21, 22, 23, 28, 29, 30, 31},
    {8, 9, 10, 11, 12, 13, 14, 15, 24, 25, 26, 27, 28, 29, 30, 31},
    {16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31}};

const int idx2[5][16] = {
    {0, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30},
    {1, 1, 5, 5, 9, 9, 13, 13, 17, 17, 21, 21, 25, 25, 29, 29},
    {3, 3, 3, 3, 11, 11, 11, 11, 19, 19, 19, 19, 27, 27, 27, 27},
    {7, 7, 7, 7, 7, 7, 7, 7, 23, 23, 23, 23, 23, 23, 23, 23},
    {15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15}};

const int idx_w[5][32] = {
    {-1, 0, 0, 1, 1, 2,  2,  3,  3,  4,  4,  5,  5,  6,  6,  7,
     7,  8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13, 14, 14, 15},
    {-1, -1, 0, 1, 1, 1, 2,  3,  3,  3,  4,  5,  5,  5,  6,  7,
     7,  7,  8, 9, 9, 9, 10, 11, 11, 11, 12, 13, 13, 13, 14, 15},
    {-1, -1, -1, -1, 0, 1, 2,  3,  3,  3,  3,  3,  4,  5,  6,  7,
     7,  7,  7,  7,  8, 9, 10, 11, 11, 11, 11, 11, 12, 13, 14, 15},
    {-1, -1, -1, -1, -1, -1, -1, -1, 0, 1, 2,  3,  4,  5,  6,  7,
     7,  7,  7,  7,  7,  7,  7,  7,  8, 9, 10, 11, 12, 13, 14, 15},
    {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13, 14, 15}};

// helper for addition
void op(LS *c1, LS *c2, const LS *g1, const LS *p1, const LS *g2, const LS *p2,
        const CK *ck) {

    LS *tmp = new_gate_bootstrapping_ciphertext(ck->params);
    bootsAND(tmp, p1, g2, ck);
    bootsOR(c1, g1, tmp, ck);
    bootsAND(c2, p1, p2, ck);

    delete_gate_bootstrapping_ciphertext(tmp);
}

// c = a + b
void par::add(LS *c, const LS *a, const LS *b, const int w, const CK *ck) {

    LS *g = new_gate_bootstrapping_ciphertext_array(w, ck->params);
    LS *p = new_gate_bootstrapping_ciphertext_array(w, ck->params);
    LS *cp = new_gate_bootstrapping_ciphertext_array(w, ck->params);

#pragma omp parallel for
    for (int i = 0; i < w * 2; ++i) {
        if (i < w) {
            bootsXOR(&p[i], &a[i], &b[i], ck);
            bootsCOPY(&cp[i], &p[i], ck);
        } else {
            bootsAND(&g[i - w], &a[i - w], &b[i - w], ck);
        }
    }

    // Sklansky's parallel prefix sum
    for (int i = 0; i < ceil(log2(w - 1)); i++) {
#pragma omp parallel for
        for (int j = 0; j <= idx_w[i][w - 1]; j++) {
            op(&g[idx1[i][j]], &p[idx1[i][j]], &g[idx1[i][j]], &p[idx1[i][j]],
               &g[idx2[i][j]], &p[idx2[i][j]], ck);
        }
    }

    bootsCOPY(&c[0], &cp[0], ck);

#pragma omp parallel for
    for (int i = 1; i < w; ++i) {
        bootsXOR(&c[i], &cp[i], &g[i - 1], ck);
    }

    delete_gate_bootstrapping_ciphertext_array(w, g);
    delete_gate_bootstrapping_ciphertext_array(w, p);
    delete_gate_bootstrapping_ciphertext_array(w, cp);
}

void par::compare(LS *c, const LS *a, const LS *b, const int w, const CK *ck) {

    LS *tmp = new_gate_bootstrapping_ciphertext(ck->params);
    LS *p = new_gate_bootstrapping_ciphertext_array(w, ck->params);
    LS *d = new_gate_bootstrapping_ciphertext_array(w, ck->params);

#pragma omp parallel for
    for (int i = 0; i < w; i++) {
        bootsXNOR(&p[i], &a[i], &b[i], ck);
        bootsCOPY(&d[i], &a[i], ck);
    }

    bootsNOT(&d[w - 1], &d[w - 1], ck);

    for (int i = 1; i < w; i *= 2) {
#pragma omp parallel for
        for (int j = 0; j < w; j += i) {
            if ((j / i) % 2 == 1) {
                bootsAND(&p[j - i], &p[j - i], &p[j], ck);
            } else {
                if (j + i < w) {
                    bootsMUX(&d[j], &p[j + i], &d[j], &d[j + i], ck);
                }
            }
        }
    }

    bootsOR(c, &d[0], &p[0], ck);

    delete_gate_bootstrapping_ciphertext(tmp);
    delete_gate_bootstrapping_ciphertext_array(w, p);
    delete_gate_bootstrapping_ciphertext_array(w, d);
}

void par::select(LS *c, const LS *s, const LS *a, const LS *b, int w,
                 const CK *ck) {
#pragma omp parallel for
    for (int i = 0; i < w; i++) {
        bootsMUX(c + i, s, a + i, b + i, ck);
    }
}

void par::select_zero(LS *c, const LS *s, const LS *a, int w, const CK *ck) {
#pragma omp parallel for
    for (int i = 0; i < w; i++) {
        bootsAND(&c[i], s, &a[i], ck);
    }
}

/*
void par::equals(LS *c, const LS *a, const LS *b,
                 const int w, const CK *ck) {

    LS *tmp = new_gate_bootstrapping_ciphertext_array(w, ck->params);

#pragma omp parallel for
    for (int i = 0; i < w; i++) {
        bootsXNOR(&tmp[i], &a[i], &b[i], ck);
    }

    for (int j = 0; (1 << j) < w; j++) {
#pragma omp parallel for
        for (int i = 0; i < w; i = i + (1 << (j + 1))) {
            if (i + (1 << j) < w) {
                bootsAND(&tmp[i], &tmp[i], &tmp[i + (1 << j)], ck);
            }
        }
    }

    bootsCOPY(r, &tmp[0], ck);

    delete_gate_bootstrapping_ciphertext_array(w, tmp);
}
*/
