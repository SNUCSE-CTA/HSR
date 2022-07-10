#include "enc.h"

int base_to_int(string a, int i) {
    if (a.length() > i) {
        if (a[i] == 'C') {
            return 1;
        }

        if (a[i] == 'G') {
            return 2;
        }

        if (a[i] == 'T') {
            return 3;
        }
    }

    return 0;
}

int alt_to_int(string a, int j) {
    int result = base_to_int(a, j * 3);
    result = (result << 2) + base_to_int(a, j * 3 + 1);
    result = (result << 2) + base_to_int(a, j * 3 + 2);
    return result;
}

void read_file(vcf *V, int line_num, string file_name) {
    ifstream infile(file_name);
    string line;
    string alt;

    for (int i = 0; i < line_num; i++) {
        if (getline(infile, line)) {
            istringstream iss(line);
            if (!(iss >> V[i].pos >> V[i].reflen >> V[i].altlen)) {
                cout << "!!!" << endl;
                break;
            } else {
                if (iss >> alt && alt[0] != '.') {
                    V[i].alt = alt;
                } else {
                    V[i].alt = "";
                }
            }
        }
    }
}

void encrypt(cvcf &X, cvcf &Y, string xfile, string yfile, int size, int n,
             int w_alt, Encryptor enc, PublicKeyPack keypack) {
    vcf x[size];
    vcf y[size];

    read_file(x, size, xfile);
    read_file(y, size, yfile);

    // int w_alt = calc_w_alt(x, y, n);
    int w_alt_size = ceil(w_alt / 6.0);

    Message x_pos(n);
    Message x_reflen(n);
    Message x_altlen(n);

    Message y_pos(n);
    Message y_reflen(n);
    Message y_altlen(n);

    vector<Message> x_alt(w_alt_size, Message(n));
    vector<Message> y_alt(w_alt_size, Message(n));

    // even slots for var_score
    for (int i = 0; i < size; i++) {

        x_pos[2 * i + 1] = {(double)x[i].pos, 0};
        x_reflen[2 * i + 1] = {(double)x[i].reflen, 0};
        x_altlen[2 * i + 1] = {(double)x[i].altlen, 0};

        y_pos[2 * i + 1] = {(double)y[i].pos, 0};
        y_reflen[2 * i + 1] = {(double)y[i].reflen, 0};
        y_altlen[2 * i + 1] = {(double)y[i].altlen, 0};

        for (int j = 0; j < w_alt_size; j++) {
            x_alt[j][2 * i + 1] = {(double)alt_to_int(x[i].alt, j), 0};
            y_alt[j][2 * i + 1] = {(double)alt_to_int(y[i].alt, j), 0};
        }
    }

    Ciphertext tmp;

    enc.encrypt(x_pos, keypack.getEncKey(), X.pos);
    enc.encrypt(x_reflen, keypack.getEncKey(), X.reflen);
    enc.encrypt(x_altlen, keypack.getEncKey(), X.altlen);

    enc.encrypt(y_pos, keypack.getEncKey(), Y.pos);
    enc.encrypt(y_reflen, keypack.getEncKey(), Y.reflen);
    enc.encrypt(y_altlen, keypack.getEncKey(), Y.altlen);

    for (int j = 0; j < w_alt_size; j++) {
        enc.encrypt(x_alt[j], keypack.getEncKey(), tmp);
        X.alt.push_back(tmp);
        enc.encrypt(y_alt[j], keypack.getEncKey(), tmp);
        Y.alt.push_back(tmp);
    }
}