
#include "complex.h"
#include <vector>
#include <iomanip>



int main(){
    cout.precision(10);
    srand(time(0));
    vector<int> log_N;
    log_N = {12, 13, 14, 15};
    vector<int> log_Q;
    log_Q = {54, 109, 219, 443}; 
    vector<int> log_P;
    log_P = {40, 40, 40, 40};
    int loop = 10;
    //instance = 0 -> complex to complex, instance = 1 -> real to real
    int instance = 1;

    for(int i = 1; i< log_N.size();i++){
        cout << "Over " << loop << " loops, the results are... \n\n" << endl;
        cout << "The parameters are: log(N) = " << log_N[i] << ", log(Q) = " << log_Q[i] << ", Delta = " << log_P[i] << "." << endl;
        //instance (0 or 1) determines whether we run complex to complex or real to real
        complex_average(log_N[i], log_Q[i], log_P[i], 1, loop, instance); 
    }
    return 0;
}