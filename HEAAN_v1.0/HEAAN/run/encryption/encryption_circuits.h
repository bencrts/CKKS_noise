/////////////////////////////////////
///the different circuits we test ///
/////////////////////////////////////

#include "../../src/HEAAN.h"
#include <NTL/ZZXFactoring.h>
#include "ntlfunctions.h"
#include "WCR_heuristics.h"

using namespace std;
using namespace NTL;

vector<double> test_encdec(long logN, long logQ, long logp)
{
    /*
	A function which determines the error due to encryption/
	decryption. 

	@param logN logarithm of the ring dimension N
	@param logQ logarithm of the (maximal) ciphertext modulus q 
	@param logp logarithm of the base used for scaling
	@return a vector of statistics about the error sizes
	*/

    // Create the context and set up the encryption scheme
    Context context(logN, logQ);
    SecretKey sk(logN);
    Scheme scheme(sk, context);
    // we always use N/2 slots
    long slots = (1 << (logN - 1));

    // Generate some random (real) data to be operated on
    double *m1 = EvaluatorUtils::randomRealArray(slots, 100);
    double *m2 = EvaluatorUtils::randomRealArray(slots, 100);

    Plaintext r = context.encode(m1, slots, logp);
    Plaintext r2 = context.encode(m2, slots, logp);

    // extract the underlying polynomial
    ZZX s = r.mx;
    ZZX s2 = r2.mx;

    // Encrypt the random arrays
    Ciphertext c1 = scheme.encrypt(m1, slots, logp, logQ);
    Ciphertext c2 = scheme.encrypt(m2, slots, logp, logQ);

    // perform some basic operations
    Ciphertext c_add = scheme.add(c1, c2);
    Ciphertext c_mult = scheme.mult(c1, c2);

    Ciphertext c_multmult = scheme.mult(c_mult, c1);
    Ciphertext c_multmultmult = scheme.mult(c_mult, c2);

    // Decrypt the ciphertexts
    Plaintext m_add = scheme.decrypt_no_decode(sk, c_add);
    Plaintext m_mult = scheme.decrypt_no_decode(sk, c_mult);

    complex<double> *m_multmult = scheme.decrypt(sk, c_multmult);
    complex<double> *m_multmultmult = scheme.decrypt(sk, c_multmultmult);

    ZZX s_add = m_add.mx;
    ZZX s_mult = m_mult.mx;
    ZZX s_adddif = s_add - s - s2;
    ZZX s_multdiff = s_mult - s * s2;

    ZZ x_z;
    ZZ y_z;
    x_z = scheme.ZZX_oo_difference(s_add, s + s2, logQ);
    double x = NTL::log(x_z);
    y_z = scheme.ZZX_oo_difference(s_mult, s * s2, logQ);
    double y = NTL::log(y_z);

    complex<double> *m3 = scheme.decode(m_add);
    complex<double> *m4 = scheme.decode(m_mult);

    // now decode and compute norms

    double size_add = 0, size_mult = 0, size_multmult = 0, size_multmultmult = 0;

    vector<double> results;

    for (int i = 0; i < slots; i++)
    {
        size_add += fabs(imag(m3[i]));
        size_mult += fabs(imag(m4[i]));
        size_multmult += fabs(imag(m_multmult[i]));
        size_multmultmult += fabs(imag(m_multmultmult[i]));
    }

    results.push_back(size_add / slots);
    results.push_back(size_mult / slots);
    results.push_back(size_multmult / slots);
    results.push_back(size_multmultmult / slots);
    
    // delete
    delete[] m1;
    delete[] m2;

    return results;
}

double basic_circuit_add(long logN, long logQ, long logp)
{
    // Create the context and set up the encryption scheme
    Context context(logN, logQ);
    SecretKey sk(logN);
    Scheme scheme(sk, context);

    //Initialise random seed
    srand(time(NULL));
    unsigned seed_1 = rand() % 10 + 1;
    unsigned seed_2 = rand() % 10 + 1;
    unsigned seed_3 = rand() % 10 + 1;

    // we always use N/2 slots
    long slots = (1 << (logN - 1));

    // bound on the oo-norm of the input vectors
    double bound = 1;

    // Generate some random (real) data to be operated on
    double *m1 = EvaluatorUtils::randomRealArray(slots, bound, seed_1);
    double *m2 = EvaluatorUtils::randomRealArray(slots, bound, seed_2);

    //Now encode
    Plaintext r1 = context.encode(m1, slots, logp);
    Plaintext r2 = context.encode(m2, slots, logp);

    // Encrypt the random arrays
    Ciphertext c1 = scheme.encrypt(m1, slots, logp, logQ);
    Ciphertext c2 = scheme.encrypt(m2, slots, logp, logQ);

    // Extract the underlying polynomial
    ZZX s1 = r1.mx;
    ZZX s2 = r2.mx;

    //Evaluate the circuit described above
    ZZX sum = s1 + s2;

    // Now do the same on the cts
    Ciphertext c_sum = scheme.add(c1, c2);

    //Now decrypt but don't decode
    Plaintext m_sum = scheme.decrypt_no_decode(sk, c_sum);

    //Extract the polynomial
    ZZX decrypted = m_sum.mx;

    //Compute the norm of their difference
    ZZ test = scheme.ZZX_oo_difference(decrypted, sum, logQ);

    double log_test = NTL::log(test);

    //Return (the log of) said difference
    ZZ two;
    two = 2;    
    // delete
    delete[] m1;
    delete[] m2;

    return log_test / NTL::log(two);
}

double basic_circuit_mult(long logN, long logQ, long logp)
{
    // Create the context and set up the encryption scheme
    Context context(logN, logQ);
    SecretKey sk(logN);
    Scheme scheme(sk, context);

    //Initialise random seed
    srand(time(NULL));
    unsigned seed_1 = rand() % 10 + 1;
    unsigned seed_2 = rand() % 10 + 1;
    unsigned seed_3 = rand() % 10 + 1;

    // we always use N/2 slots
    long slots = (1 << (logN - 1));

    // bound on the oo-norm of the input vectors
    double bound = 1;

    // Generate some random (real) data to be operated on
    double *m1 = EvaluatorUtils::randomRealArray(slots, bound, seed_1);
    double *m2 = EvaluatorUtils::randomRealArray(slots, bound, seed_2);

    //Now encode
    Plaintext r1 = context.encode(m1, slots, logp);
    Plaintext r2 = context.encode(m2, slots, logp);

    // Encrypt the random arrays
    Ciphertext c1 = scheme.encrypt(m1, slots, logp, logQ);
    Ciphertext c2 = scheme.encrypt(m2, slots, logp, logQ);

    // Extract the underlying polynomial
    ZZX s1 = r1.mx;
    ZZX s2 = r2.mx;

    //Evaluate the circuit described above
    ZZX prod = s1 * s2;

    // Now do the same on the cts
    Ciphertext c_mult = scheme.mult(c1, c2);

    //Now decrypt but don't decode
    Plaintext m_mult = scheme.decrypt_no_decode(sk, c_mult);

    //Extract the polynomial
    ZZX decrypted = m_mult.mx;

    //Compute the norm of their difference
    ZZ test = scheme.ZZX_oo_difference(decrypted, prod, logQ);

    double log_test = NTL::log(test);

    //Return (the log of) the difference
    ZZ two;
    two = 2;    

    // delete
    delete[] m1;
    delete[] m2;

    return log_test / NTL::log(two);
}

double encoded_circuit(long logN, long logQ, long logp)
{
    /*A function that generates three random numbers
	m1, m2, m3. It then encodes those values into ri.
	On the one side, it computes prod = (r1 + r2)*r3.
	On the other, it encrypts the values into ci. Then
	computes decrypted = decrypt_no_decode((c1 + c2)*c3).
	We then compute the difference prod - dcrypted.
	This is to measure the size of the LWE error. 

	@param logN logarithm of the ring dimension N
	@param logQ logarithm of the (maximal) ciphertext modulus q 
	@param logp logarithm of the base used for scaling
	@return a vector of statistics about the error sizes
	*/

    // Create the context and set up the encryption scheme
    Context context(logN, logQ);
    SecretKey sk(logN);
    Scheme scheme(sk, context);

    //Initialise random seed
    srand(time(NULL));
    unsigned seed_1 = rand() % 10 + 1;
    unsigned seed_2 = rand() % 10 + 1;
    unsigned seed_3 = rand() % 10 + 1;

    // we always use N/2 slots
    long slots = (1 << (logN - 1));

    // bound on the oo-norm of the input vectors
    double bound = 1;

    // Generate some random (real) data to be operated on
    double *m1 = EvaluatorUtils::randomRealArray(slots, bound, seed_1);
    double *m2 = EvaluatorUtils::randomRealArray(slots, bound, seed_2);
    double *m3 = EvaluatorUtils::randomRealArray(slots, bound, seed_3);

    //Now encode
    Plaintext r1 = context.encode(m1, slots, logp);
    Plaintext r2 = context.encode(m2, slots, logp);
    Plaintext r3 = context.encode(m3, slots, logp);

    ZZ x = polynomial_max_coeff(r1.mx);

    //cout << x << endl;

    // Encrypt the random arrays
    Ciphertext c1 = scheme.encrypt(m1, slots, logp, logQ);
    Ciphertext c2 = scheme.encrypt(m2, slots, logp, logQ);
    Ciphertext c3 = scheme.encrypt(m3, slots, logp, logQ);

    // Extract the underlying polynomial
    ZZX s1 = r1.mx;
    ZZX s2 = r2.mx;
    ZZX s3 = r3.mx;

    //Evaluate the circuit described above
    ZZX sum = s1 + s2;
    ZZX prod = sum * s3;
    ZZX plain_rs;
    Ring2Utils::rightShift(plain_rs, prod, logp, pow(2, logN));

    // Now do the same on the cts
    Ciphertext c_add = scheme.add(c1, c2);
    Ciphertext c_mult = scheme.mult(c_add, c3);

    //Rescale
    scheme.reScaleByAndEqual(c_mult, logp);

    //Now decrypt but don't decode
    Plaintext m_decrypted = scheme.decrypt_no_decode(sk, c_mult);

    //Extract the polynomial
    ZZX decrypted = m_decrypted.mx/;

    //Compute the norm of their difference
    ZZ test = scheme.ZZX_oo_difference(decrypted, prod, logQ);

    double log_test = NTL::log(test);

    //Return the difference
    ZZ two;
    two = 2;    
    // delete
    delete[] m1;
    delete[] m2;
    delete[] m3;

    return log_test / NTL::log(two);
}


