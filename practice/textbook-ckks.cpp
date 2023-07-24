#include <vector>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "../HEAAN/HEAAN/src/HEAAN.h"

using namespace std;
using namespace NTL;

// a debugging function for checking on randomness. Remove from final.
void printVector(complex<double> *z, int count)
{
    cout << "(";
    for (int i = 0; i < count; i++)
    {
        cout << z[i].real() << " + " << z[i].imag() << "i ,";
    }
    cout << " ... )\n";
}

/**
 * @brief adds the array first `slots` entries of the array z2 references to the first `slots` entries of the array z1 references.
 * Will throw out of bounds if size mismatch or `slots` too large.
 *
 * @tparam T type which permits +=
 * @param z1 pointer to the array to add to
 * @param z2 pointer to the array to be added
 * @param slots number of entries to add
 */
template <typename T>

void addInplace(T *z1, T *z2, long slots)
{
    for (int i = 0; i < slots; i++)
    {
        z1[i] += z2[i];
    }
}

/**
 * @brief multiplies the array first `slots` entries of the array z1 references by the first `slots` entries of the array z2 references.
 * Will throw out of bounds if size mismatch or `slots` too large.
 *
 * @tparam T type which permits *=
 * @param z1 pointer to the array to add to
 * @param z2 pointer to the array to be added
 * @param slots number of entries to add
 */
template <typename T>
void multiplyInplace(T *z1, T *z2, long slots)
{
    for (int i = 0; i < slots; i++)
    {
        z1[i] *= z2[i];
    }
}

// this is from helib/src/norms.cpp
ZZ largestCoeff(const ZZX &f)
{
    ZZ mx = ZZ::zero();
    for (long i = 0; i <= deg(f); i++)
    {
        if (mx < abs(coeff(f, i)))
            mx = abs(coeff(f, i));
    }
    return mx;
}

/**
 * @brief calculate the infinity norm of both the real and complex difference.
 *
 * @param z1 pointer to the first array of complex numbers
 * @param z2 pointer to the second array of complex numbers
 * @param slots number of entries
 * @return array<double, 2> norm of real difference, norm of imaginary difference
 */
array<double, 2> difference(const complex<double> *z1, const complex<double> *z2, long slots)
{
    array<double, 2> diff{{0, 0}};
    for (int i = 0; i < slots; i++)
    {
        double r_diff = abs(z1[i].real() - z2[i].real());
        double c_diff = abs(z1[i] - z2[i]);
        diff[0] = (diff[0] > r_diff) ? diff[0] : r_diff;
        diff[1] = (diff[1] > c_diff) ? diff[1] : c_diff;
    }
    return diff;
}

ZZ ringDifference(complex<double> *z1, complex<double> *z2, Context context, long slots, long logDelta)
{
    // this doesn't account for encoding error if one of the z's does encode without rounding
    ZZX m1 = context.encode(z1, slots, logDelta);
    ZZX m2 = context.encode(z2, slots, logDelta);
    return largestCoeff(m1 - m2);
}

/**
 * @brief run the circuit `loop` times with the given parameters and data with real and imaginary components sampled from [0, input_bound].
 * Report on the maximum and average noise in the infinity norm in three spaces: `plaintext`, `real`, and `complex`.
 *
 * @param log_N Polynomial degree.
 * @param log_Q Top level plaintext modulus. In this library, also gives the size of the auxiliary keyswitch modulus.
 * @param log_Delta Precision parameter.
 * @param input_bound Bound on real and imaginary components of the randomly sampled input data.
 * @param loop How many trials.
 */
void noiseExperiment(long logN, long logQ, long logDelta, long h, int input_bound, int loop)
{
    // entires are laid out real_average, real_max, complex_average, complex_max
    array<double, 4> add_noises{{0, 0, 0, 0}};
    array<double, 4> mult_noises{{0, 0, 0, 0}};
    array<ZZ, 2> ring_add_noises{{NTL::ZZ::zero(), NTL::ZZ::zero()}};
    array<ZZ, 2> ring_mult_noises{{NTL::ZZ::zero(), NTL::ZZ::zero()}};

    for (int j = 0; j < loop; j++)
    {
        // initialise the scheme.
        SecretKey sk(logN, h);
        Context context(logN, logQ);
        Scheme scheme(sk, context);
        long slots = (1 << (logN - 1));

        // generate three random complex vectors.
        complex<double> *z1 = EvaluatorUtils::randomComplexArray(slots, input_bound);
        complex<double> *z2 = EvaluatorUtils::randomComplexArray(slots, input_bound);
        complex<double> *z3 = EvaluatorUtils::randomComplexArray(slots, input_bound);

        // Now encrypt
        Ciphertext c1 = scheme.encrypt(z1, slots, logDelta, logQ);
        Ciphertext c2 = scheme.encrypt(z2, slots, logDelta, logQ);
        Ciphertext c3 = scheme.encrypt(z3, slots, logDelta, logQ);

        // plaintext addition: z1 += z2
        addInplace(z1, z2, slots);

        // ciphertext addition: c1 += c2
        scheme.addAndEqual(c1, c2);

        // plaintext multiplication: z3 *= z1
        multiplyInplace(z3, z1, slots);

        // ciphertext multiplication: c3 *= c1;
        scheme.multAndEqual(c3, c1);
        scheme.reScaleByAndEqual(c3, logDelta);

        // first decrypt
        complex<double> *z_add = scheme.decrypt(sk, c1);
        complex<double> *z_mult = scheme.decrypt(sk, c3);

        // find both real and complex differences
        array<double, 2> add_noise = difference(z_add, z1, slots);
        array<double, 2> mult_noise = difference(z_mult, z3, slots);

        // update msg space noises: arranged real_average, real_max, complex_average, complex_max
        add_noises[0] += add_noise[0];
        add_noises[1] = (add_noise[0] > add_noises[1]) ? add_noise[0] : add_noises[1];
        add_noises[2] += add_noise[1];
        add_noises[3] = (add_noise[1] > add_noises[3]) ? add_noise[1] : add_noises[3];

        mult_noises[0] += mult_noise[0];
        mult_noises[1] = (mult_noise[0] > mult_noises[1]) ? mult_noise[0] : mult_noises[1];
        mult_noises[2] += mult_noise[1];
        mult_noises[3] = (mult_noise[1] > mult_noises[3]) ? mult_noise[1] : mult_noises[3];

        // find ring space differences
        ZZ ring_add_noise = ringDifference(z_add, z1, context, slots, logDelta);
        ZZ ring_mult_noise = ringDifference(z_mult, z3, context, slots, logDelta);

        // and update ring space noises: arranged average, max
        ring_add_noises[0] += ring_add_noise;
        ring_add_noises[1] = (ring_add_noise > ring_add_noises[1]) ? ring_add_noise : ring_add_noises[1];
        ring_mult_noises[0] += ring_mult_noise;
        ring_mult_noises[1] = (ring_mult_noise > ring_mult_noises[1]) ? ring_mult_noise : ring_mult_noises[1];

        delete[] z1, z2, z3, z_add, z_mult;
    }

    // find msg space averages on entries 0 and 2 by dividing by loop
    add_noises[0] /= loop;
    add_noises[2] /= loop;
    mult_noises[0] /= loop;
    mult_noises[2] /= loop;


    // legible prints
    cout << "The parameters are: log(N) = " << logN << ", log(Q) = " << logQ << ", log(Delta) = " << logDelta << ", h = " << h << ".\n";
    cout << "over " << loop << " trials, the results are:" << endl;
    cout << "REAL ADDITION, average = " << add_noises[0] << ", maximum = " << add_noises[1] << endl;
    cout << "REAL MULTIPLICATION, average = " << mult_noises[0] << ", maximum = " << mult_noises[1] << endl;
    cout << "COMPLEX ADDITION, average = " << add_noises[2] << ", maximum = " << add_noises[3] << endl;
    cout << "COMPLEX MULTIPLICATION, average = " << mult_noises[2] << ", maximum = " << mult_noises[3] << endl;
    cout << "RING ADDITION, average = " << conv<RR>(ring_add_noises[0]) / loop << ", maximum = " << ring_add_noises[1] << endl;
    cout << "RING MULTIPLICATION, average = " << conv<RR>(ring_mult_noises[0]) / loop << ", maximum = " << ring_mult_noises[1] << endl;
}

int main()
{
    srand(1222);                                  // seed for message generation
    NTL::SetSeed(NTL::conv<NTL::ZZ>((long)1222)); // seed for encryption randomness

    array<int, 3> log_N{{13, 14, 15}};
    array<int, 3> log_Q{{109, 219, 443}};
    array<int, 3> log_Delta{{40, 40, 40}};
    array<int, 3> h{{64, 64, 64}};
    array<int, 3> input_bound{{1, 1, 1}};

    int loop = 10;

    for (int i = 0; i < log_N.size(); i++)
    {
        noiseExperiment(log_N[i], log_Q[i], log_Delta[i], h[i], input_bound[i], loop);
    }

    return 0;
}
