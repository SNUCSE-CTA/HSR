#include "enc.h"
#include <cmath>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace std;

void toArray(int *result, int a, int len) {
    for (int i = 0; i < len; i++) {
        result[i] = (a >> i) & 1;
    }
}

void alt_to_bit(int *a, string alt_str, int w_alt) {
    for (int j = 0; j < w_alt; j++) {
        a[j] = 0;
    }

    for (int j = 0; j < alt_str.length() && j < w_alt / 2; j++) {
        if (alt_str[j] == 'C') {
            a[2 * j] = 0;
            a[2 * j + 1] = 1;
        } else if (alt_str[j] == 'G') {
            a[2 * j] = 1;
            a[2 * j + 1] = 0;
        } else if (alt_str[j] == 'T') {
            a[2 * j] = 1;
            a[2 * j + 1] = 1;
        }
    }
}

void read_vcf(vcf *V, string file_name, int line_num, int w_alt) {
    ifstream infile(file_name);
    string line;
    string alt_str;
    for (int i = 0; i < line_num; i++) {
        if (getline(infile, line)) {
            istringstream iss(line);
            if (!(iss >> V[i].pos >> V[i].reflen >> V[i].altlen)) {
                cout << "!!!" << endl;
                break;
            } else {
                // V[i].alt = new int[w_alt];
                if (iss >> alt_str) {
                    alt_to_bit(V[i].alt, alt_str, w_alt);
                }
            }
        }
    }
}

void encrypt(cvcf *S, string file_name, const int size, const int w,
             const int w_alt, const TFheGateBootstrappingSecretKeySet *key,
             const TFheGateBootstrappingParameterSet *params) {

    int tmp[w];
    vcf a[size];

    read_vcf(a, file_name, size, w_alt);

    for (int i = 0; i < size; i++) {

        S[i].pos = new_gate_bootstrapping_ciphertext_array(w, params);
        S[i].reflen = new_gate_bootstrapping_ciphertext_array(w, params);
        S[i].altlen = new_gate_bootstrapping_ciphertext_array(w, params);
        S[i].alt = new_gate_bootstrapping_ciphertext_array(w_alt, params);

        // encrypt POS
        toArray(tmp, a[i].pos, w);
        for (int j = 0; j < w; j++) {
            bootsSymEncrypt(&S[i].pos[j], tmp[j], key);
        }

        // encrypt REFLEN
        toArray(tmp, a[i].reflen, w);
        for (int j = 0; j < w; j++) {
            bootsSymEncrypt(&S[i].reflen[j], tmp[j], key);
        }

        // encrypt ALTLEN
        toArray(tmp, a[i].altlen, w);
        for (int j = 0; j < w; j++) {
            bootsSymEncrypt(&S[i].altlen[j], tmp[j], key);
        }

        // encrypt ALT
        for (int j = 0; j < w_alt; j++) {
            bootsSymEncrypt(&S[i].alt[j], a[i].alt[j], key);
        }
    }
}

void encrypt_to_file(string file_name, const int size, const int w,
                     const int w_alt,
                     const TFheGateBootstrappingSecretKeySet *key,
                     const TFheGateBootstrappingParameterSet *params) {

    int tmp[w];
    vcf a[size];
    size_t lastindex = file_name.find_last_of(".");
    string out_file_name = file_name.substr(0, lastindex);
    out_file_name += ".enc";

    read_vcf(a, file_name, size, w_alt);

    FILE *cloud_data = fopen(out_file_name.c_str(), "wb");

    LS *S = new_gate_bootstrapping_ciphertext_array(3 * w + w_alt, params);

    for (int i = 0; i < size; i++) {
        // encrypt POS
        for (int j = 0; j < w; j++) {
            toArray(tmp, a[i].pos, w);
            bootsSymEncrypt(&S[j], tmp[j], key);
            export_gate_bootstrapping_ciphertext_toFile(cloud_data, &S[j],
                                                        params);
        }
        // encrypt REFLEN
        for (int j = 0; j < w; j++) {
            toArray(tmp, a[i].reflen, w);
            bootsSymEncrypt(&S[j + w], tmp[j], key);
            export_gate_bootstrapping_ciphertext_toFile(cloud_data, &S[j + w],
                                                        params);
        }
        // encrypt ALTLEN
        for (int j = 0; j < w; j++) {
            toArray(tmp, a[i].altlen, w);
            bootsSymEncrypt(&S[j + 2 * w], tmp[j], key);
            export_gate_bootstrapping_ciphertext_toFile(cloud_data,
                                                        &S[j + 2 * w], params);
        }

        // encrypt ALT
        for (int j = 0; j < w_alt; j++) {
            bootsSymEncrypt(&S[j + 3 * w], a[i].alt[j], key);
            export_gate_bootstrapping_ciphertext_toFile(cloud_data,
                                                        &S[j + 3 * w], params);
        }
    }

    fclose(cloud_data);
    delete_gate_bootstrapping_ciphertext_array(3 * w + w_alt, S);
}