#include "enc.h"
#include "eval.h"
#include "omp.h"
#include <chrono>
#include <iostream>

using namespace std;

int decrypt_val(LS *a, int w, TFheGateBootstrappingSecretKeySet *seckey) {
    int mask = -1;

    int val = 0;
    int ai = 0;
    for (int k = 0; k < w; k++) {
        ai = bootsSymDecrypt(&a[k], seckey);
        if (ai == 1) {
            val += 1 << k;
        }
    }

    if (ai == 1) {
        val = (val | (mask << w));
    }

    return val;
}

int main(int argc, char *argv[]) {

    // Default params
    int size = 5;   // |X| = |Y|
    int R_len = 27; // |R|
    int w_int = 9;  // w_int
    int w_alt = 24; // w_alt

    // Scoring scheme
    int s_m = 5;
    int s_s = -3;
    int g_o = -9;
    int g_e = -1;

    string file_x = "../dataset/paper_x.dat"; // list X
    string file_y = "../dataset/paper_y.dat"; // list Y

    // Parse params
    if (argc == 6) {
        file_x = argv[1];
        file_y = argv[2];
        size = atoi(argv[3]);
        R_len = atoi(argv[4]);
        w_int = atoi(argv[5]);
        // w_alt = atoi(argv[6]);
    }

    cout << "Input" << endl << "-----------------" << endl;
    cout << "X: " << file_x << endl;
    cout << "Y: " << file_y << endl;
    cout << "|X|=|Y|: " << size << endl;
    cout << "|R|: " << R_len << endl;
    cout << "w_int: " << w_int << endl;
    cout << "w_alt: " << w_alt << endl;
    cout << "Scpring scheme: {" << s_m << "/" << s_s << "," << g_o << "," << g_e
         << "}" << endl;
    cout << "-----------------" << endl;

    // TFHE key and parameter geneation
    const int minimum_lambda = 128;

    TFheGateBootstrappingParameterSet *params =
        new_default_gate_bootstrapping_parameters(minimum_lambda);

    uint32_t seed[] = {314, 1592, 1093};
    tfhe_random_generator_setSeed(seed, 3);
    TFheGateBootstrappingSecretKeySet *key =
        new_random_gate_bootstrapping_secret_keyset(params);

    FILE *cloud_key = fopen("cloud.key", "wb");
    export_tfheGateBootstrappingCloudKeySet_toFile(cloud_key, &key->cloud);
    fclose(cloud_key);

    FILE *cloud_key2 = fopen("cloud.key", "rb");
    CK *ck = new_tfheGateBootstrappingCloudKeySet_fromFile(cloud_key2);
    fclose(cloud_key2);

    // Encryption
    cvcf X[size];
    cvcf Y[size];

    cout << "Encrypting input ..." << endl;

    encrypt(X, file_x, size, w_int, w_alt, key, params);
    encrypt(Y, file_y, size, w_int, w_alt, key, params);

    cout << "Encryption done." << endl;

    // Initialize variables
    LS *V[size];
    LS *M[size + 1];
    LS *S[size + 1];
    LS *E[size + 1];

    for (int i = 0; i < size; i++) {
        V[i] = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
        M[i] = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
        S[i] = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
        E[i] = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    }

    M[size] = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    S[size] = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    E[size] = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);

    LS *m = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *s = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);
    LS *e = new_gate_bootstrapping_ciphertext_array(w_int, ck->params);

    cout << "Computing scores of regions ..." << endl;

    auto start = chrono::steady_clock::now();
    score(V, M, S, E, X, Y, size, R_len, s_m, s_s, g_o, g_e, w_int, w_alt, ck);
    auto end = chrono::steady_clock::now();

    auto diff = end - start;
    int sec = chrono::duration<double, milli>(diff).count() / 1000;
    cout << "Elapsed time (s): " << sec << endl;

    /*
    //Print V and M
    for (int i = 0; i < size; i++) {
        cout << decrypt_val(V[i], w_int, key) << " ";
    }
    std::cout << endl;

    for (int i = 0; i < size + 1; i++) {
        cout << decrypt_val(M[i], w_int, key) << " ";
    }
    std::cout << endl;
    */

    cout << "Finding highly similar region ..." << endl;
    auto start2 = chrono::steady_clock::now();
    find(m, s, e, V, M, S, E, size, w_int, ck);
    auto end2 = chrono::steady_clock::now();

    auto diff2 = end2 - start2;
    int sec2 = chrono::duration<double, milli>(diff2).count() / 1000;
    cout << "Elapsed time (s): " << sec2 << endl;

    cout << "Decrypting output ..." << endl;
    int m_dec = decrypt_val(m, w_int, key);
    int s_dec = decrypt_val(s, w_int, key);
    int e_dec = decrypt_val(e, w_int, key);
    cout << "Decryption done. " << endl;

    cout << "Output" << endl << "-----------------" << endl;
    cout << "Score: " << m_dec << endl;
    cout << "Start pos: " << s_dec << endl;
    cout << "End pos: " << e_dec << endl;

    // Free memory
    delete_gate_bootstrapping_secret_keyset(key);
    delete_gate_bootstrapping_cloud_keyset(ck);
    delete_gate_bootstrapping_parameters(params);

    for (int i = 0; i < size; i++) {
        delete_gate_bootstrapping_ciphertext_array(w_int, V[i]);
        delete_gate_bootstrapping_ciphertext_array(w_int, M[i]);
        delete_gate_bootstrapping_ciphertext_array(w_int, S[i]);
        delete_gate_bootstrapping_ciphertext_array(w_int, E[i]);
    }

    delete_gate_bootstrapping_ciphertext_array(w_int, M[size]);
    delete_gate_bootstrapping_ciphertext_array(w_int, S[size]);
    delete_gate_bootstrapping_ciphertext_array(w_int, E[size]);

    delete_gate_bootstrapping_ciphertext_array(w_int, m);
    delete_gate_bootstrapping_ciphertext_array(w_int, s);
    delete_gate_bootstrapping_ciphertext_array(w_int, e);

    return 0;
}