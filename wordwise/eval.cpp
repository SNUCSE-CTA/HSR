#include "eval.h"

double len_scale = 100.0;
double score_scale = 1000.0;

// Gap penalty function (affine)
void W(Ciphertext &r, Ciphertext l, int g_o, int g_e, Ciphertext one,
       HomEvaluator eval, PublicKeyPack keypack) {

    Ciphertext t1, t2;
    eval.sub(l, one, t1);
    eval.mult(t1, g_e, t2);
    eval.add(t2, g_o, r);
    // cout << "W() " << endl;
    // cout << "r " << t2.getModulusBits() << " -> ";
    MathUtils::approxDiscreteEqualZero(eval, keypack, l, t1);
    // cout << t1.getModulusBits() << endl;
    eval.sub(one, t1, t2);
    eval.mult(t2, r, keypack.getMultKey(), r);
}

// HSR_W-VarScore
void var_score(Ciphertext &V, Ciphertext &lref, cvcf X, cvcf Y, int n, int s_m,
               int s_s, int g_o, int g_e, int w_alt, int d, Ciphertext one,
               HomEvaluator eval, PublicKeyPack keypack) {

    Ciphertext m, g, gx, gy, t1, t2;

    // Eq. (7)
    eval.mult(X.reflen, 1.0 / len_scale, t1);
    eval.mult(Y.reflen, 1.0 / len_scale, t2);
    // cout << "Lref " << t1.getModulusBits() << " -> ";
    MathUtils::approxMinMax(eval, keypack, t1, t2, t1, lref, d);
    // cout << lref.getModulusBits() << endl;
    eval.mult(lref, len_scale, lref);

    // Eq. (8) and (9) (Begin)
    eval.sub(X.altlen, X.reflen, t1);
    eval.sub(Y.altlen, Y.reflen, t2);

    // Eq. (12) and (13)
    eval.add(lref, t1, gx);
    eval.add(lref, t2, gy);

    // Eq. (8) and (9) (Cont), t1 = min, t2 = max
    eval.mult(t1, 1.0 / len_scale, t1);
    eval.mult(t2, 1.0 / len_scale, t2);
    // cout << "lmin " << t1.getModulusBits() << " -> ";
    MathUtils::approxMinMax(eval, keypack, t1, t2, t1, t2, d);
    // cout << t1.getModulusBits() << endl;
    eval.mult(t1, len_scale, t1);
    eval.mult(t2, len_scale, t2);

    // Eq. (6)
    eval.sub(t2, t1, g);
    // cout << "g: " << g.getModulusBits() << endl;
    // Eq. (5)
    eval.add(lref, t1, m);

    // Checking x.alt == y.alt
    MathUtils::approxDiscreteEqual(eval, keypack, X.altlen, Y.altlen, t1);
    // cout << "e: " << t1.getModulusBits() << " -> ";

    int w_alt_len = ceil(w_alt / 6.0);

    for (int j = 0; j < w_alt_len; j++) {
        MathUtils::approxDiscreteEqual(eval, keypack, X.alt[j], Y.alt[j], t2);
        eval.mult(t1, t2, keypack.getMultKey(), t1);
    }
    // cout << t1.getModulusBits() << endl;

    Ciphertext v1, v2;

    // Eq. (10)
    W(v1, g, g_o, g_e, one, eval, keypack);
    eval.mult(m, s_m, t2);
    eval.mult(t2, t1, keypack.getMultKey(), t2);
    eval.add(v1, t2, v1);

    eval.sub(one, t1, t1);
    eval.mult(m, s_s, t2);
    eval.mult(t2, t1, keypack.getMultKey(), t2);
    eval.add(v1, t2, v1);

    // Eq. (14)
    W(t1, gx, g_o, g_e, one, eval, keypack);
    W(t2, gy, g_o, g_e, one, eval, keypack);
    eval.add(t1, t2, v2);

    // cout << "v1 " << v1.getModulusBits() << endl;
    // cout << "v2 " << v2.getModulusBits() << endl;

    // Eq. (15)
    eval.mult(v1, 1.0 / score_scale, t1);
    eval.mult(v2, 1.0 / score_scale, t2);
    MathUtils::approxMinMax(eval, keypack, t1, t2, t1, t2, d);
    eval.mult(t2, score_scale, V);

    // cout << "V " << V.getModulusBits() << endl;
}

// HSR_W-Score
void score(Ciphertext &V, Ciphertext &M, cvcf X, cvcf Y, int n, int R_len,
           int s_m, int s_s, int g_o, int g_e, int w_alt, int d, Ciphertext one,
           HomEvaluator eval, PublicKeyPack keypack) {

    Ciphertext S, E, lref, t;

    var_score(V, lref, X, Y, n, s_m, s_s, g_o, g_e, w_alt, d, one, eval,
              keypack);

    // int N = X.reflen.getNumberOfSlots();
    Message P(n);
    P[0] = -1;
    P[2 * n] = R_len + 1;

    eval.add(lref, X.pos, t);
    eval.rightRotate(t, 1, keypack.getRightRotKey(1), S);
    eval.leftRotate(X.pos, 1, keypack.getLeftRotKey(1), t);
    eval.add(t, P, E);
    eval.sub(E, S, t);
    eval.mult(t, s_m, M);
}

// HSR_W-Find
void find(Ciphertext &score, Ciphertext V, Ciphertext M, int size, int n, int d,
          HomEvaluator eval, PublicKeyPack keypack) {

    int log_n = ceil(log2(2 * size + 1)) + 1;

    Message mask(n);
    Message minf(n);

    for (int i = 0; i < 2 * size + 2; i++) {
        mask[i] = {1.0, 0};
    }

    for (int i = 2 * size + 2; i < n; i++) {
        minf[i] = {-1.0, 0.0};
    }

    Ciphertext t1, psmax, psmax_scaled;

    eval.add(V, M, t1);
    eval.mult(t1, mask, t1);
    eval.rightRotate(t1, 1, keypack.getRightRotKey(1), psmax);

    // cout << "MV " << psmax.getModulusBits() << endl;
    // print(psmax, N, dec, secret_key);

    for (int i = 0; i < log_n - 1; i++) {
        eval.rightRotate(psmax, 1 << i, keypack.getRightRotKey(1 << i), t1);
        // cout << "RR " << i << endl;
        // print(t1, N, dec, secret_key);
        eval.add(psmax, t1, psmax);
        // cout << "ADDED " << i << endl;
        // print(psmax, N, dec, secret_key);
    }

    // cout << "PSUM: " << psmax.getModulusBits() << endl;
    // print(psmax, N, dec, secret_key);

    eval.mult(psmax, mask, t1); // t2 stores scaled down psmax
    eval.mult(t1, 1.0 / score_scale, psmax_scaled);
    eval.add(psmax_scaled, minf, psmax);

    for (int i = 0; i < log_n - 1; i++) {
        // cout << i << endl;
        // cout << psmax.getModulusBits() << " -> " << endl;
        eval.leftRotate(psmax, 1 << i, keypack.getLeftRotKey(1 << i), t1);
        MathUtils::approxMinMax(eval, keypack, psmax, t1, t1, psmax, d);
        // cout << psmax.getModulusBits() << endl;
    }

    // cout << "smax: " << endl;
    // printOne(psmax, 2 * n + 1, dec, secret_key);

    eval.sub(psmax, psmax_scaled, t1);
    eval.mult(t1, mask, psmax);

    // cout << "sub: " << endl;
    // printOne(psmax, 0, dec, secret_key);

    // find maximum of cmax
    for (int i = 0; i < log_n - 1; i++) {
        // cout << i << endl;
        // cout << psmax.getModulusBits() << " -> " << endl;
        eval.leftRotate(psmax, 1 << i, keypack.getLeftRotKey(1 << i), t1);
        MathUtils::approxMinMax(eval, keypack, psmax, t1, t1, psmax, d);
        // cout << psmax.getModulusBits() << endl;
    }

    // cout << "max: " << endl;
    // printOne(psmax, 0, dec, secret_key);

    // scale down
    eval.mult(psmax, score_scale, score);
}
