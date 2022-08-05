#include "noise-final.h"
#include <vector>

int main()
{

    srand(time(0));

    vector<int> log_N;
    log_N = {12, 13, 14, 15};
    vector<int> log_Q;
    log_Q = {54, 109, 219, 443}; // prev values: 109, 218, 438, 886
    vector<int> log_P;
    log_P = {40, 40, 40, 40};
    ZZ plaintext_bound = conv<ZZ>("1099511627776");
    int loop = 10;

    int instance = 1;

    if (instance == 0)
    {
        cout << "We are running noise_final_fixed_vector: \n\n";
    }
    else if (instance == 1)
    {
        cout << "We are running noise_final_random_plaintext: \n\n";
    }

    for (int i = 0; i < log_N.size(); i++)
    {
        cout << "Over " << loop << " loops, the results are... \n\n"
             << endl;
        cout << "The parameters are: log(N) = " << log_N[i] << ", log(Q) = " << log_Q[i] << ", Delta = " << log_P[i] << ".\n"
             << endl;

        // replace in noise-final.h with appropriate function
        // this will determine if we run: random reals, random PTs or fixed poly
        SecretKey sk(log_N[i]);
        noise_final_average(log_N[i], log_Q[i], log_P[i], plaintext_bound, loop, sk, instance);
    }
    return 0;
}
