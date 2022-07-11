#include "enc.h"
#include "eval.h"

using namespace std;
using namespace heaan;

vector<double> decrypt_val(Ciphertext ctxt, int len, Decryptor dec,
                           SecretKey sk) {
    Message dmsg;
    dec.decrypt(ctxt, sk, dmsg);

    vector<double> result(len);
    for (long i = 0; i < len; i++) {
        result.push_back(dmsg[i].real());
    }

    return result;
}

int main(int argc, char *argv[]) {

    int thread_num = thread::hardware_concurrency();
    SetNumThreads(thread_num);
    srand(time(NULL));

    // Scoring scheme
    double s_m = 5;
    double s_s = -3;
    int g_o = -9;
    int g_e = -1;

    // Scheme params
    int l = 15;
    int q = 28;
    int d = 17;
    int w_alt = 24;

    int size = 18;
    string file_x = "../dataset/paper_x.dat"; // list X
    string file_y = "../dataset/paper_y.dat"; // list Y
    int R_len = 123;

    if (argc == 8) {
        file_x = argv[1];
        file_y = argv[2];
        size = atoi(argv[3]);
        R_len = atoi(argv[4]);
        l = atoi(argv[5]);
        q = atoi(argv[6]);
        d = atoi(argv[7]);
    }

    // score_scale = (double)(rgen << 3);
    // score_scale = 10;
    // while (score_scale < s_m * R_len) {
    //    score_scale *= 10;
    //}

    int log_n = ceil(log2(2 * size + 1)) + 1;
    int n = 1 << log_n;

    cout << "Input" << endl << "-----------------" << endl;
    cout << "X: " << file_x << endl;
    cout << "Y: " << file_y << endl;
    cout << "|X|=|Y|: " << size << endl;
    cout << "|R|: " << R_len << endl;
    cout << "l: " << l << endl;
    cout << "q: " << q << endl;
    cout << "d: " << d << endl;
    cout << "n: " << n << endl;
    // cout << "Score scale: " << score_scale << endl;

    cout << "Generating keys and parameters ..." << endl;
    // HEAAN key and parameter generation
    Parameters params(l, q);
    Context context(params);
    context.makeBootstrappable(log_n);

    string key_dir_path = "./keys";
    mkdir(key_dir_path.data(), 0775);

    SecretKey secret_key(context);
    PublicKeyPack keypack(context, secret_key, key_dir_path);

    Encryptor enc(context);
    Decryptor dec(context);
    HomEvaluator eval(context);

    cout << "Generation done." << endl;
    cout << "HEAAN Params" << endl;
    cout << "----------------" << endl;
    cout << "N: " << params.getLogDegree() << endl;
    // cout << "q: " << params.getQuantizeBits() << endl;
    cout << "Modulus bits: " << params.getModulusBits() << endl;
    cout << "----------------" << endl;

    cvcf X, Y;
    cout << "Encrypting input ..." << endl;
    encrypt(X, Y, file_x, file_y, size, n, w_alt, enc, keypack);
    cout << "Encryption done." << endl;

    // vector<double> ds = decrypt_val(X.pos, n, dec, secret_key);
    // for (double b : ds)
    //     cout << b << " ";
    // cout << endl;

    // Ciphertext of (1, 1, ..., 1)
    Message m_one(n);
    for (int i = 0; i < 2 * size + 1; i++) {
        m_one[i] = {1.0, 0};
    }

    Ciphertext one;
    enc.encrypt(m_one, keypack.getEncKey(), one);

    cout << "Computing scores of regions ..." << endl;

    Ciphertext V, M, sub;
    auto start = chrono::steady_clock::now();
    score(V, M, X, Y, n, R_len, s_m, s_s, g_o, g_e, w_alt, d, one, eval,
          keypack);
    auto end = chrono::steady_clock::now();

    auto diff = end - start;
    int sec = chrono::duration<double, milli>(diff).count() / 1000;
    cout << "Elapsed time (s): " << sec << endl;

    cout << "Finding highly similar region ..." << endl;
    auto start2 = chrono::steady_clock::now();
    find(sub, V, M, size, n, d, eval, keypack);
    auto end2 = chrono::steady_clock::now();

    auto diff2 = end2 - start2;
    int sec2 = chrono::duration<double, milli>(diff2).count() / 1000;
    cout << "Elapsed time (s): " << sec2 << endl;

    cout << "Decrypting output ..." << endl;
    vector<double> score_dec = decrypt_val(sub, n, dec, secret_key);
    cout << "Decryption done. " << endl;

    cout << "Output" << endl << "-----------------" << endl;
    cout << "Score: " << score_dec[0] << endl;

    return 0;
}
