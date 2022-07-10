#include "seq_circuit.h"

// c = a + b
void seq::add(LS *c, const LS *a, const LS *b, const int w, const CK *ck) {

    LS *carry = new_gate_bootstrapping_ciphertext_array(2, ck->params);
    LS *temp = new_gate_bootstrapping_ciphertext_array(3, ck->params);

    bootsCONSTANT(carry, 0, ck);

    for (int i = 0; i < w; ++i) {
        bootsXOR(temp, a + i, b + i, ck);
        bootsXOR(c + i, temp, carry, ck);
        bootsAND(temp + 1, a + i, b + i, ck);
        bootsAND(temp + 2, carry, temp, ck);
        bootsXOR(carry + 1, temp + 1, temp + 2, ck);
        bootsCOPY(carry, carry + 1, ck);
    }

    // bootsCOPY(&sum[w], &carry[0], ck); //last carry = overflow

    delete_gate_bootstrapping_ciphertext_array(3, temp);
    delete_gate_bootstrapping_ciphertext_array(2, carry);
}

// c = (a >= b) 1:0
void seq::compare(LS *c, const LS *a, const LS *b, const int w, const CK *ck) {

    LS *tmp = new_gate_bootstrapping_ciphertext_array(2, ck->params);

    bootsCONSTANT(tmp, 1, ck);

    // comparison
    for (int i = 0; i < w - 1; i++) {
        bootsXNOR(tmp + 1, a + i, b + i, ck);
        bootsMUX(tmp, tmp + 1, tmp, a + i, ck);
    }

    // sign bit comparison
    bootsXNOR(tmp + 1, a + w - 1, b + w - 1, ck);
    bootsNOT(c, a + w - 1, ck);
    bootsMUX(c, tmp + 1, tmp, c, ck);

    delete_gate_bootstrapping_ciphertext_array(2, tmp);
}

void seq::enc_const(LS *c, const int a, int w, const CK *ck) {
    for (int i = 0; i < w; i++) {
        bootsCONSTANT(&c[i], (a >> i) & 1, ck);
    }
}

// c = (a==b)? 1: 0;
void seq::equals(LS *c, const LS *a, const LS *b, const int w, const CK *ck) {

    LS *tmp = new_gate_bootstrapping_ciphertext_array(w, ck->params);

    for (int i = 0; i < w; i++) {
        bootsXNOR(&tmp[i], a + i, b + i, ck);
    }

    bootsCOPY(c, &tmp[0], ck);

    for (int i = 1; i < w; i++) {
        bootsAND(c, c, &tmp[i], ck);
    }

    delete_gate_bootstrapping_ciphertext_array(w, tmp);
}

// c = (a==0)? 1: 0;
void seq::equals_zero(LS *c, const LS *a, const int w, const CK *ck) {

    bootsCONSTANT(c, 0, ck);

    for (int i = 0; i < w; ++i) {
        bootsOR(c, c, &a[i], ck);
    }

    bootsNOT(c, c, ck);
}

// c = a - b
void seq::subtract(LS *c, const LS *a, const LS *b, const int w, const CK *ck) {

    LS *carry = new_gate_bootstrapping_ciphertext_array(2, ck->params);

    bootsCONSTANT(&carry[0], 1, ck);

    LS *temp = new_gate_bootstrapping_ciphertext_array(3, ck->params);

    for (int i = 0; i < w; i++) {
        bootsXNOR(temp, a + i, b + i, ck);
        bootsXOR(c + i, temp, carry, ck);
        bootsANDYN(temp + 1, a + i, b + i, ck);
        bootsAND(temp + 2, carry, temp, ck);
        bootsXOR(carry + 1, temp + 1, temp + 2, ck);
        bootsCOPY(carry, carry + 1, ck);
    }

    // bootsCOPY(&r[w], &carry[0], ck); //last carry = overflow

    delete_gate_bootstrapping_ciphertext_array(3, temp);
    delete_gate_bootstrapping_ciphertext_array(2, carry);
}

// c = a - 1
void seq::subtract_one(LS *c, const LS *a, const int w, const CK *ck) {

    LS *carry = new_gate_bootstrapping_ciphertext(ck->params);

    bootsCONSTANT(carry, 1, ck);

    for (int i = 0; i < w; i++) {
        bootsXOR(c + i, a + i, carry, ck);
        bootsANDNY(carry, a + i, carry, ck);
    }

    delete_gate_bootstrapping_ciphertext(carry);
}

void seq::select(LS *c, const LS *s, const LS *a, const LS *b, int w,
                 const CK *ck) {
    for (int i = 0; i < w; i++) {
        bootsMUX(c + i, s, a + i, b + i, ck);
    }
}

void seq::select_zero(LS *c, const LS *s, const LS *a, const int w,
                      const CK *ck) {
    for (int i = 0; i < w; i++) {
        bootsAND(&c[i], s, &a[i], ck);
    }
}

// c = a + 1
void add_one(LS *c, const LS *a, const int w, const CK *ck) {

    LS *carry = new_gate_bootstrapping_ciphertext(ck->params);

    bootsCONSTANT(carry, 1, ck);

    for (int i = 0; i < w; ++i) {
        bootsXOR(c + i, a + i, carry, ck);
        bootsAND(carry, a + i, carry, ck);
    }

    delete_gate_bootstrapping_ciphertext(carry);
}

// c = a << bitnum
// void shift_left(LS *c, int bitnum, const LS *a, int w, const CK *ck) {
//     seq::copy(c + bitnum, a, w - bitnum, ck);
//     seq::fill_zero(c, bitnum, ck);
// }

// c = a * t
void seq::multiply_const(LS *c, const LS *a, const int t, const int w,
                         const CK *ck) {
    LS *tmp = new_gate_bootstrapping_ciphertext_array(w, ck->params);
    LS *ca = new_gate_bootstrapping_ciphertext_array(w, ck->params);

    copy(ca, a, w, ck);
    fill_zero(c, w, ck);

    int abst = abs(t);
    int log = ceil(log2(abst + 1));

    for (int i = 0; i < log; i++) {
        if (abst % 2 == 1) {
            add(tmp, c, ca, w, ck);
            copy(c, tmp, w, ck);
        }

        // shift_left(tmp, 1, ca, w, ck);
        copy(tmp + 1, ca, w - 1, ck);
        fill_zero(tmp, 1, ck);

        copy(ca, tmp, w, ck);
        abst = abst / 2;
    }

    if (t < 0) {
        // negate(c, c, w, ck);
        for (int i = 0; i < w; i++) {
            bootsNOT(tmp + i, c + i, ck);
        }
        add_one(c, tmp, w, ck);
    }

    delete_gate_bootstrapping_ciphertext_array(w, tmp);
}

void seq::copy(LS *a, const LS *b, const int w, const CK *ck) {
    for (int i = 0; i < w; i++) {
        bootsCOPY(a + i, b + i, ck);
    }
}

void seq::fill_zero(LS *a, int w, const CK *ck) {
    for (int i = 0; i < w; i++) {
        bootsCONSTANT(a + i, 0, ck);
    }
}
