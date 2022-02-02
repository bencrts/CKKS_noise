#include "NR_circuit.h"
#include "../../src/HEAAN.h"
#include <math.h>
#include <NTL/ZZXFactoring.h>
#include <stdlib.h>
using namespace std;
using namespace NTL;

int main()
{
    int loop = 10;
    //NR stuff
    double a = 1;
    double b = 5;
    int iternum = 5;
    //HEAAN stuff
    long logN = 12;
    //sorry the P is too confusing for me
    long logDelta = 25;
    NR_circuit(a,b,iternum,loop,logN,logDelta);

    return 0;
}
