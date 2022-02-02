#include "NRfunctions.h"
#include "functions.h"
#include "../../src/HEAAN.h"

using namespace std;
using namespace NTL;

void NR_circuit(double a,double b,int iternum,int loop,long logN, long logDelta){
    long logQ = 10 + 2 * logDelta * (1 + iternum);
    string filename = "NR|"+to_string(round(a))+"|"+to_string(round(b))+"|"+to_string(iternum)+"|"+to_string(loop)+"|"+to_string(logN)+"|"+to_string(logDelta);
    //containers for error
    dVec max;
    max.assign(iternum + 1, -RAND_MAX);
    dVec avg;
    avg.assign(iternum + 1, 0.0);

    Context context(logN, logQ);
    long slots = (1 << (logN - 1));
    cVec cresult;
    double *rresult;
    double maxdiff;

    for (int j = 0; j < loop; j++)
    {
        SecretKey sk(logN);
        Scheme scheme(sk, context);

        //generate some data
        double *z = random_array(slots, a, b);
        print_array(slots, z, 10);
        print_array(slots, fz(slots, z), 10);

        Ciphertext zc = scheme.encrypt(z, slots, logDelta, logQ);
        //we need our zc each iteration, so make a copy which will update
        Ciphertext fzc = zc;

        //iteration 0
        scheme.multByConstAndEqual(fzc, T1(a, b), logDelta);
        scheme.reScaleByAndEqual(fzc, logDelta);
        scheme.addConstAndEqual(fzc, T0(a, b), logDelta);

        cresult = scheme.decrypt(sk, fzc);
        rresult = real_part(slots, cresult);

        //find correct bits on iter 0
        maxdiff = max_log_diff(slots, fz(slots, z), rresult);
        cout << maxdiff << "\n";
        update_max_avg(max, avg, 0, maxdiff);
        Ciphertext fztemp;
        for (int i = 0; i < iternum; i++)
        {
            fztemp = scheme.add(fzc, fzc);

            scheme.squareAndEqual(fzc);
            scheme.reScaleByAndEqual(fzc, logDelta);

            scheme.multAndEqual(fzc, zc);
            scheme.reScaleByAndEqual(fzc, logDelta);

            scheme.negateAndEqual(fzc);
            scheme.addAndEqual(fzc, fztemp);

            cresult = scheme.decrypt(sk, fzc);
            rresult = real_part(slots, cresult);

            double maxdiff = max_log_diff(slots, fz(slots, z), rresult);
            cout << maxdiff << "\n";
            update_max_avg(max, avg, i + 1, maxdiff);
        }
        delete[] z;
    }
    scale_vec(avg,loop);
    print_vec_toFile(filename, max);
    print_vec_toFile(filename, avg);
}