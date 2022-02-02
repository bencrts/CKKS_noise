// final noise, hopefully
// copies from new-stuff.h and editing

#include "../../src/HEAAN.h"
#include <math.h>
#include <stdlib.h>
#include "../../src/Ring2Utils.h"
#include "../encryption/ntlfunctions.h"
//#include "new-stuff.h"

using namespace std;
using namespace NTL;

Plaintext attained_polynomial(long logN, long B, long logP, long logQ)
{
    /* Generates a plaintext for ring dimension 2**logN, with
    coefficients all equal to B

    @param logN the log of the ring dimension
    @param B the value of the plaintext coefficients
    */

    ZZX f_out;

    for (int i = 0; i < pow(2, logN); i++)
    {
        SetCoeff(f_out, i, B);
    }

    Plaintext p; //= Plaintext(f_out,logP,logQ,logN -1,false);
    p.mx = f_out;
    p.logp = logP;
    p.logq = logQ;
    p.slots = logN - 1;
    p.isComplex = false;

    return p;
}

vector<double> complex_to_complex(long logN, long logQ, long logp, double bound = 1.0)
{
    SecretKey sk(logN);
    Context context(logN, logQ);
    Scheme scheme(sk, context);
    long slots = (1 << (logN - 1));

    // Generate some random (real) data to be operated on
    complex<double> *m1 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double> *m2 = EvaluatorUtils::randomComplexArray(slots, bound);
    complex<double> *m3 = EvaluatorUtils::randomComplexArray(slots, bound);

    // Now encrypt
    Ciphertext c1 = scheme.encrypt(m1, slots, logp, logQ);
    Ciphertext c2 = scheme.encrypt(m2, slots, logp, logQ);
    Ciphertext c3 = scheme.encrypt(m3, slots, logp, logQ);
    vector<double> precision(2);

    // Now let's operate on raw data
    complex<double> *add_raw = complex_vector_add(m1, m2, slots);
    complex<double> *mult_raw = complex_vector_mult(m3, add_raw, slots);

    // Now add, and mult
    // keyswitch done modulo qQ, actual "multiplication" done modulo q (the "current" ciphertext modulus)
    Ciphertext c_add = scheme.add(c1, c2);
    Ciphertext c_mult_no_rescale = scheme.mult(c3, c_add);

    // rescale 
    Ciphertext c_mult = scheme.reScaleBy(c_mult_no_rescale, logp);

    // decrypting
    complex<double> *decrypted_add = scheme.decrypt(sk, c_add);
    complex<double> *decrypted_mult = scheme.decrypt(sk, c_mult);

    // now let's get the precision loss
    double precision_loss_difference_add = max_difference(add_raw, decrypted_add, slots);
    double precision_loss_difference_mult = max_difference(mult_raw, decrypted_mult, slots);

    precision[0] = precision_loss_difference_add;
    precision[1] = precision_loss_difference_mult;

    return precision;
}

vector<double> real_to_complex(long logN, long logQ, long logp, double bound = 1.0)
{
    SecretKey sk(logN);
    Context context(logN, logQ);
    Scheme scheme(sk, context);
    long slots = (1 << (logN - 1));

    // Generate some random (real) data to be operated on
    double *m1 = EvaluatorUtils::randomRealArray(slots, bound);
    double *m2 = EvaluatorUtils::randomRealArray(slots, bound);
    double *m3 = EvaluatorUtils::randomRealArray(slots, bound);

    // convert these to complex arrays
    complex<double> *cm1 = real_to_complex(m1, slots);
    complex<double> *cm2 = real_to_complex(m2, slots);
    complex<double> *cm3 = real_to_complex(m3, slots);

    // Now encrypt
    Ciphertext c1 = scheme.encrypt(m1, slots, logp, logQ);
    Ciphertext c2 = scheme.encrypt(m2, slots, logp, logQ);
    Ciphertext c3 = scheme.encrypt(m3, slots, logp, logQ);

    delete[] m1, m2, m3;
    vector<double> precision(2);

    // Now let's operate on raw data
    complex<double> *add_raw = complex_vector_add(cm1, cm2, slots);
    complex<double> *mult_raw = complex_vector_mult(cm3, add_raw, slots);

    // Now add, and mult
    // keyswitch done modulo qQ, actual "multiplication" done modulo q (the "current" ciphertext modulus)
    Ciphertext c_add = scheme.add(c1, c2);
    Ciphertext c_mult_no_rescale = scheme.mult(c3, c_add);

    // rescale
    Ciphertext c_mult = scheme.reScaleBy(c_mult_no_rescale, logp);

    // decrypt
    complex<double> *decrypted_add = scheme.decrypt(sk, c_add);
    complex<double> *decrypted_mult = scheme.decrypt(sk, c_mult);

    // now let's get the precision loss
    double precision_loss_real_difference_add = max_real_difference(add_raw, decrypted_add, slots);
    double precision_loss_real_difference_mult = max_real_difference(mult_raw, decrypted_mult, slots);

    precision[0] = precision_loss_real_difference_add;
    precision[1] = precision_loss_real_difference_mult;

    // delete
    delete[] cm1, cm2, cm3, add_raw, mult_raw, decrypted_add, decrypted_mult;
    return precision;
}

void complex_average(long logN, long logQ, long logp, long bound, int loop, int instance)
{
    vector<double> results(2);
    vector<double> max_values(2);
    for (int i = 0; i < loop; i++)
    {
        vector<double> noise;
        switch (instance)
        {
        case 0:
            noise = complex_to_complex(logN, logQ, logp, bound);
            break;
        case 1:
            noise = real_to_complex(logN, logQ, logp, bound);
            break;
        }

        results[0] += noise[0];
        results[1] += noise[1];

        for (int j = 0; j < 2; j++)
        {
            if (noise[j] > max_values[j])
            {
                max_values[j] = noise[j];
            }
        }
    }
    for (int j = 0; j < 2; j++)
    {
        results[j] = results[j] / loop;
    }

    switch (instance)
    {
    case 0:
        cout << "We are running random complex numbers. \n\n";
        break;
    case 1:
        cout << "We are running random real numbers. \n\n";
        break;
    }
    cout << "The average noise for addition is: " << results[0] << ".\n";
    cout << "The max noise for addition is: " << max_values[0] << ".\n\n";
    cout << "The average noise for multiplication is: " << results[1] << ".\n";
    cout << "The max noise for multiplication is: " << max_values[1] << ".\n\n";
    cout << "in bits,\n";
    cout << "The average noise for addition is: " << log2(results[0]) << " bits\n";
    cout << "The max noise for addition is: " << log2(max_values[0]) << " bits\n\n";
    cout << "The average noise for multiplication is: " << log2(results[1]) << " bits\n";
    cout << "The max noise for multiplication is: " << log2(max_values[1]) << " bits.\n\n";
    cout << "Latex +: " << "\n";
    cout << logN << " & " << logQ << " & " << log2(results[0]) << " & " << log2(max_values[0]) << " &  &  & " << "\n";
    cout << "Latex x: " << endl;
    cout << logN << " & " << logQ << " & " << log2(results[1]) << " & " << log2(max_values[1]) << " &  &  & " << "\n";
}
