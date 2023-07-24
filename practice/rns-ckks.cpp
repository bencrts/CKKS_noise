#include <string>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <ctime>

#include "../FullRNS-HEAAN/src/Ciphertext.h"
#include "../FullRNS-HEAAN/src/Common.h"
#include "../FullRNS-HEAAN/src/Context.h"
#include "../FullRNS-HEAAN/src/EvaluatorUtils.h"
#include "../FullRNS-HEAAN/src/Key.h"
#include "../FullRNS-HEAAN/src/Numb.h"
#include "../FullRNS-HEAAN/src/Plaintext.h"
#include "../FullRNS-HEAAN/src/Scheme.h"
#include "../FullRNS-HEAAN/src/SchemeAlgo.h"
#include "../FullRNS-HEAAN/src/SecretKey.h"
#include "../FullRNS-HEAAN/src/StringUtils.h"
#include "../FullRNS-HEAAN/src/TimeUtils.h"

using namespace std;
complex<double> *convert_real_to_complex(double *array, long slots)
{
    complex<double> *complex_array = new complex<double>[slots];
    for (int i = 0; i < slots; i++)
    {
        complex_array[i].real(array[i]);
        complex_array[i].imag(0);
    }
    return complex_array;
}

double* convert_complex_to_real(complex<double>* complex_array, long slots){

	double* real_array = new double[slots];
	for (int i = 0; i < slots; i++)
    {
        real_array[i] = real(complex_array[i]);
    }

	return real_array;
}

double max_difference(complex<double> *array1, complex<double> *array2, long slots)
{
    double max = 0;
    complex<double> diff;
    double size;
    for (int i = 0; i < slots; i++)
    {
        diff = array1[i] - array2[i];
        size = abs(diff);
        if (size > max)
        {
            max = size;
        }
    }
    return max;
}

double max_real_difference(complex<double> *array1, complex<double> *array2, long slots)
{
    double max = 0;
    complex<double> diff;
    double size;
    for (int i = 0; i < slots; i++)
    {
        // grab real part of the complex arrays
        diff = array1[i].real() - array2[i].real();
        size = abs(diff);
        if (size > max)
        {
            max = size;
        }
    }
    return max;
}

vector<double> complex_to_complex_and_real(long logN, long L, long logp, long logSlots, long k){

    srand(time(0));
    // set up the scheme, leaving all parameters to be determined as function inputs
    Context context(logN, logp, L, k);
    SecretKey secretKey(context);
    Scheme scheme(secretKey, context);
    long slots = (1 << logSlots);
    double bound = 1.0;

    // generate random complex data with a bound of 1.0 (as before)
    complex<double>* mvec1 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvec2 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvec3 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double>* mvecAdd = new complex<double>[slots];
    complex<double>* mvecMult = new complex<double>[slots];

    // associated plaintext operations to compare against
    for(long i = 0; i < slots; i++) {
        mvecAdd[i] = mvec1[i] + mvec2[i];
        mvecMult[i] = mvecAdd[i] * mvec3[i];
    }

    // encrypt
    Ciphertext cipher1 = scheme.encrypt(mvec1, slots, L);
    Ciphertext cipher2 = scheme.encrypt(mvec2, slots, L);
    Ciphertext cipher3 = scheme.encrypt(mvec3, slots, L);

    // operate (addition, multplication and rescale)
    Ciphertext addCipher = scheme.add(cipher1, cipher2);
    Ciphertext multCipher = scheme.mult(addCipher, cipher3);
    scheme.reScaleByAndEqual(multCipher, 1);

    // decrypt
    complex<double>* dvecAdd = scheme.decrypt(secretKey, addCipher);
    complex<double>* dvecMult = scheme.decrypt(secretKey, multCipher);

    // compute precision loss
	double real_precision_loss_difference_add = max_real_difference(dvecAdd, mvecAdd, slots);
	double real_precision_loss_difference_mult = max_real_difference(dvecMult, mvecMult, slots);
    double complex_precision_loss_difference_add = max_difference(dvecAdd, mvecAdd, slots);
    double complex_precision_loss_difference_mult = max_difference(dvecMult, mvecMult, slots);

    vector<double> precision(4);
    precision[0] = real_precision_loss_difference_add;
    precision[1] = real_precision_loss_difference_mult;
    precision[2] = complex_precision_loss_difference_add;
    precision[3] = complex_precision_loss_difference_mult;
    return precision;
}


void real_and_complex_average(long logN, long L, long logp, long logSlots, long k, int loop)
{
    vector<double> results(4);
    vector<double> max_values(4);
    for (int i = 0; i < loop; i++)
    {
        vector<double> noise;

		noise = complex_to_complex_and_real(logN, L, logp, logSlots, k);

    results[0] += noise[0];
    results[1] += noise[1];
		results[2] += noise[2];
		results[3] += noise[3];

        for (int j = 0; j < 4; j++)
        {
            if (noise[j] > max_values[j])
            {
                max_values[j] = noise[j];
            }
        }
    }
    for (int j = 0; j < 4; j++)
    {
        results[j] = results[j] / loop;
    }

    cout << "The parameters are: log(N) = " << logN << ", log(Q) = " << 60 + logp * (L-1) << ", log(Delta) = " << logp << ".\n";
    cout << "over " << loop << " trials, the results are:" << endl;
    cout << "REAL ADDITION, average = " << results[0] << ", maximum = " << max_values[0] << endl;
    cout << "REAL MULTIPLICATION, average = " << results[1] << ", maximum = " << max_values[1] << endl;
    cout << "COMPLEX ADDITION, average = " << results[2] << ", maximum = " << max_values[2] << endl;
    cout << "COMPLEX MULTIPLICATION, average = " << results[3] << ", maximum = " << max_values[3] << endl;
    cout << "\n\n";

}


int main(){

  int loop = 10;
	real_and_complex_average(12, 2, 40, 11, 3, loop);
	real_and_complex_average(13, 2, 40, 12, 3, loop);
	real_and_complex_average(14, 5, 40, 13, 6, loop);
	real_and_complex_average(15, 10, 40, 14, 11, loop);

}
