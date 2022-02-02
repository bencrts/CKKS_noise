/////////////////////////////////////
///NTL specific utility functions ///
/////////////////////////////////////

#include <NTL/ZZXFactoring.h>
#include <vector>
#include <iostream>
#include <random>

#include "../../src/HEAAN.h"

using namespace std;
using namespace NTL;

void print_results(vector<double> results)
{
    /*
	A function which prints the elements of a results vector
	@param results a vector of results
	*/
    cout << "printing result..." << endl;
    for (int i = 0; i < results.size(); i++)
    {
        cout << results[i] << endl;
    }
}

void print_complex_array(complex<double> *array, int count)
{
    /*
	A function which prints the elements of a results vector
	@param results a vector of results
	*/
    //cout << "printing result..." << endl;
    for (int i = 0; i < count - 1; i++)
    {
        cout << array[i].real() << "+" << array[i].imag() << "i,";
    }
    cout << array[count - 1].real() << "+" << array[count - 1].imag() << "i...\n\n";
}

Plaintext generate_random_plaintext(long logN, ZZ B, long logP, long logQ)
{

    /* Generates a random plaintext for ring dimension 2**logN, with
	coefficients between -B and B

	@param logN the log of the ring dimension
	@param B the bound on the plaintext coefficients

	EXAMPLE:

	ZZX g = generate_random_plaintext(2, 100);
	cout << "The polynomial is " << g << endl;

	OUTPUT:
	The polynomial is [24 0 88 -59]
	*/
    srand(time(0));

    ZZX f_out;
    ZZ x;

    for (int i = 0; i < pow(2, logN); i++)
    {
        x = RandomBnd(B);
        // this samples integers between [0, B], let's add a 50% chance this is -ve
        if (rand() % 2 == 0)
        {
            x = (-1) * x;
        }
        SetCoeff(f_out, i, x);
    }
    Plaintext p; //= Plaintext(f_out,logP,logQ,logN -1,false);
    p.mx = f_out;
    p.logp = logP;
    p.logq = logQ;
    p.slots = logN - 1;
    p.isComplex = false;
    return p;
}

complex<double> *generate_random_plaintext_v2(long logN, ZZ B, long logp, long logQ)
{
    //so first make the random poly
    srand(time(0));

    ZZX f_out;
    ZZ x;

    for (int i = 0; i < pow(2, logN); i++)
    {
        x = RandomBnd(B);
        // this samples integers between [0, B], let's add a 50% chance this is -ve
        if (rand() % 2 == 0)
        {
            x = (-1) * x;
        }
        SetCoeff(f_out, i, x);
    }

    // now map back to complex!
    long slots = 1 << (logN - 1);
    complex<double> *vec = new complex<double>[slots];
    Context context(logN, 1);
    ZZ tmp;
    for (long i = 0, idx = 0; i < slots; ++i, ++idx)
    {

        // real part of res[i]
        tmp = f_out[idx];
        vec[i].real(EvaluatorUtils::scaleDownToReal(tmp, logp));

        // imaginary part of res[i]
        tmp = f_out[idx + context.Nh];
        vec[i].imag(EvaluatorUtils::scaleDownToReal(tmp, logp));
    }

    context.fftSpecial(vec, slots);
    return vec;
}

ZZ polynomial_max_coeff(ZZX f)
{
    /*
	Get the maximal coefficient of an input polynomial f, i.e. if 
	f = a_0 + a_1x + ... + a_{n-1}x^{n-1} we return the value
	a = max_i {a_i}

	@param f an NTL polynomial

	EXAMPLE:

	// f = 1 + 2x - 3x^2
	ZZX f; 
	SetCoeff(f, 0, 1);
	SetCoeff(f, 1, 2);
	SetCoeff(f, 2, -3);
	ZZ max = polynomial_max_coeff(f);
	cout << max << endl;

	OUTPUT:
	3
	*/

    ZZ max;
    max = 0;

    for (int i = 0; i < deg(f) + 1; i++)
    {
        // get the vaue a_i
        ZZ a_i = NTL::coeff(f, i);

        // make any negative coeffs positive
        if (a_i < 0)
        {
            a_i = (-1) * a_i;
        }

        // cycle through all coeffs, keeping track of the max
        if (a_i > max)
        {
            max = a_i;
        }
    }

    return max;
}

double poly_oo_difference(ZZX poly1, ZZX poly2, long logQ)
{
    ZZX diff = poly1 - poly2;
    int d = deg(diff);
    double max = 0;
    for (int i = 0; i < d + 1; i++)
    {

        if (abs(diff[i]) > max)
            max = conv<double>(abs(diff[i]));
    }
    return max;
}

void print_poly(ZZX m, int count)
{
    cout << m[0] << "+";
    for (int i = 1; i < count; i++)
    {
        cout << m[i] << "X^" << i << "+";
    }
    cout << "...\n\n";
}

complex<double> *constant_poly_vector(long logN)
{
    // want to return a vector of complex numbers that encodes to a constant polynomial
    // this shouldn't depend on logQ or logDelta, so try and force myself to write it
    // without them
    long slots = 1 << (logN - 1);
    complex<double> *vec = new complex<double>[slots]();
    for (int i = 0; i < slots; i++)
    {
        vec[i].real(1);
        vec[i].imag(1);
    }
    Context context(logN, 1);
    context.fftSpecial(vec, slots);
    //print_poly(mx, 11);
    return vec;
}

double max_difference(complex<double> *array1, complex<double> *array2, long slots)
{
    // we're taking two arrays of complex numbers and finding the infinity norm of their difference --
    // unlike max_real_difference, we are looking at the difference as a complex number & then taking the absolute value
    // code cannot check the arrays are both of the same size
    double max = 0;
    complex<double> diff;
    double size;
    for (int i = 0; i < slots; i++)
    {
        //take difference on ith slot
        diff = array1[i] - array2[i];
        //take norm: but this is squared
        size = norm(diff);
        //so squareroot
        size = sqrt(size);
        //and compare
        if (size > max)
        {
            max = size;
        }
    }
    return max;
}

double max_real_difference(complex<double> *array1, complex<double> *array2, long slots)
{
    // this time we're only looking for the maximum REAL difference i.e. max_j Re(x_j - y_j)
    // again no check on the arrays being the same length
    double max = 0;
    double rdiff;
    for (int j = 0; j < slots; j++)
    {
        //take the difference of real part
        rdiff = array1[j].real() - array2[j].real();
        //absolute value
        rdiff = abs(rdiff);
        if (rdiff > max)
        {
            max = rdiff;
        }
    }
    return max;
}

complex<double> *complex_vector_add(complex<double> *array1, complex<double> *array2, long slots)
{
    // we're going to add the two 'vectors' represented by array1 and array2 and return a new one
    // slots is the length of the vectors. No check that this is correct or that vectors match in size
    complex<double> *sum = new complex<double>[slots]();
    for (int i = 0; i < slots; i++)
    {
        sum[i] = array1[i] + array2[i];
    }
    return sum;
}

complex<double> *complex_vector_mult(complex<double> *array1, complex<double> *array2, long slots)
{
    // component wise product of two complex vectors. Slots is length of vectors
    complex<double> *product = new complex<double>[slots]();
    for (int j = 0; j < slots; j++)
    {
        product[j] = array1[j] * array2[j];
    }
    return product;
}

complex<double> *real_to_complex(double *array, long slots)
{
    // we're going to take a vector of reals and convert to a complex array so we can use our complex functions
    complex<double> *complex_array = new complex<double>[slots];
    for (int i = 0; i < slots; i++)
    {
        complex_array[i].real(array[i]);
        complex_array[i].imag(0);
    }
    return complex_array;
}