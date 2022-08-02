// final noise, hopefully
// copies from new-stuff.h and editing

#include "../../src/HEAAN.h"
#include <math.h>
#include <stdlib.h>
#include "../../src/Ring2Utils.h"
#include "ntlfunctions.h"

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

vector<ZZ> noise_final_random_reals(long logN, long logQ, long logp, long bound, SecretKey sk)
{
  // cout << "######## PARAMETERS: logN = " << logN << " logQ = " << logQ << " logDelta = " << logp << " bound = " << bound << "\n\n";
  Context context(logN, logQ);
  // SecretKey sk;
  Scheme scheme(sk, context);
  long slots = (1 << (logN - 1));

  unsigned seed_1 = rand() % 10 + 1;
  unsigned seed_2 = rand() % 10 + 1;
  unsigned seed_3 = rand() % 10 + 1;

  // Generate some random (real) data to be operated on
  double *m1 = EvaluatorUtils::randomRealArray(slots, bound, seed_1);
  double *m2 = EvaluatorUtils::randomRealArray(slots, bound, seed_2);
  double *m3 = EvaluatorUtils::randomRealArray(slots, bound, seed_3);

  // Now encode
  Plaintext r1 = context.encode(m1, slots, logp);
  Plaintext r2 = context.encode(m2, slots, logp);
  Plaintext r3 = context.encode(m3, slots, logp);

  // Now encrypt
  // I know the subscripts are a bit weird. But I wanted
  // to keep c3 as in new-stuff.h to avoid potential errors
  Ciphertext c0 = scheme.encrypt(m1, slots, logp, logQ);
  Ciphertext c1 = scheme.encrypt(m2, slots, logp, logQ);
  Ciphertext c2 = scheme.encrypt(m3, slots, logp, logQ);

  // Trying the vector
  vector<ZZ> noises(8);

  // Now add, and mult
  // logQ watch: keyswitch done modulo qQ, actual "multiplication" done modulo q (the "current" ciphertext modulus)

  Ciphertext c_add = scheme.add(c0, c1);
  Ciphertext c3 = scheme.mult(c2, c_add);
  // logQ watch: no reduction here

  // change: I'm going to do this modulo Q, the ciphertext modulus, to match the multiplication within scheme above
  ZZ q = context.qpowvec[c1.logq];
  ZZX add;  // add of the plaintexts
  ZZX prod; // prod of the plaintexts
  Ring2Utils::add(add, r1.mx, r2.mx, q, context.N);
  Ring2Utils::mult(prod, r3.mx, add, q, context.N);

  // printing add; no need for rescale, this shouldn't need anything else
  // cout << "printing r0 + r1 mod Q:\n\n";
  // print_poly(add, 12);
  Plaintext addp = add;

  // printing prod: this hasn't been rescaled yet
  // cout << "printing r1*r2 mod Q:\n\n";
  // print_poly(prod, 12);
  Plaintext prodp = prod;

  // running a decryption on c3, which is c_add*c2 without a rescale
  // also running a decryption of c_add = c_0 + c_1

  Plaintext decrypted_add = scheme.decrypt_no_decode(sk, c_add);
  Plaintext decrypted = scheme.decrypt_no_decode(sk, c3);

  // cout << "Printing the decryption of c0+c1:\n\n";
  // print_poly(decrypted_add.mx, 12);
  // cout << "Printing the decryption of c_add*c2:\n\n";
  // print_poly(decrypted.mx, 12);

  // balance the addition (plaintext)
  ZZX addb = scheme.balance_polynomial(add, decrypted_add.logq);
  // balance ciphertext add balanced (ciphertext)
  ZZX caddb = scheme.balance_polynomial(decrypted_add.mx, decrypted_add.logq);

  // balance the prod (plaintext)
  ZZX prodb = scheme.balance_polynomial(prod, decrypted.logq);
  // balance ciphertext prod balanced (ciphertext)
  ZZX cprodb = scheme.balance_polynomial(decrypted.mx, decrypted.logq);

  // now manually get the noise. Do this a variety of ways and see what shakes out
  // I think all of these are the same for small enough noise but potentially will be different once
  // the noise wraps around q/2....

  // using fi for add noise
  // ei for mult noise

  // first: just use the scheme's in built function, being careful to use the correct logq
  ZZ f1 = scheme.ZZX_oo_difference(add, decrypted_add.mx, decrypted_add.logq);
  // second: get the max coefficient of the difference
  ZZX diff_add = add - decrypted_add.mx;
  ZZ f2 = polynomial_max_coeff(diff_add);
  // third: balance the difference and get the max coefficient
  ZZX diffb_add = scheme.balance_polynomial(diff_add, decrypted_add.logq);
  ZZ f3 = polynomial_max_coeff(diffb_add);
  // fourth: take the difference mod q instead
  ZZX diffq_add;
  q = context.qpowvec[decrypted_add.logq];
  Ring2Utils::sub(diffq_add, add, decrypted_add.mx, q, context.N);
  ZZ f4 = polynomial_max_coeff(diffq_add);
  /*cout << "\n######## ADD NOISE ########\n\n";
  cout << "WITHOUT BALANCING:\n";
  cout << "noise f1 = " << f1 << "\n";
  cout << "noise f2 = " << f2 << "\n";
  cout << "noise f3 = " << f3 << "\n";
  cout << "noise f4 = " << f4 << "\n\n"; */

  // first: just use the scheme's in built function, being careful to use the correct logq
  f1 = scheme.ZZX_oo_difference(addb, caddb, decrypted_add.logq);
  // second: get the max coefficient of the difference
  diff_add = addb - caddb;
  f2 = polynomial_max_coeff(diff_add);
  // third: balance the difference and get the max coefficient
  diffb_add = scheme.balance_polynomial(diff_add, decrypted_add.logq);
  f3 = polynomial_max_coeff(diffb_add);
  // fourth: take the difference mod q instead
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq_add, addb, caddb, q, context.N);
  f4 = polynomial_max_coeff(diffq_add);
  /*cout << "ADD WITH BALANCING:\n";
  cout << "noise f1 = " << f1 << "\n";
  cout << "noise f2 = " << f2 << "\n";
  cout << "noise f3 = " << f3 << "\n";
  cout << "noise f4 = " << f4 << "\n\n"; */

  noises[0] = f1;
  noises[1] = f2;
  noises[2] = f3;
  noises[3] = f4;

  // first: just use the scheme's in built function, being careful to use the correct logq
  ZZ e1 = scheme.ZZX_oo_difference(prod, decrypted.mx, decrypted.logq);
  // second: get the max coefficient of the difference
  ZZX diff = prod - decrypted.mx;
  ZZ e2 = polynomial_max_coeff(diff);
  // third: balance the difference and get the max coefficient
  ZZX diffb = scheme.balance_polynomial(diff, decrypted.logq);
  ZZ e3 = polynomial_max_coeff(diffb);
  // fourth: take the difference mod q instead
  ZZX diffq;
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq, prod, decrypted.mx, q, context.N);
  ZZ e4 = polynomial_max_coeff(diffq);
  /*cout << "\n######## MULT NOISE BEFORE RESCALE ########\n\n";
  cout << "WITHOUT BALANCING:\n";
  cout << "noise e1 = " << e1 << "\n";
  cout << "noise e2 = " << e2 << "\n";
  cout << "noise e3 = " << e3 << "\n";
  cout << "noise e4 = " << e4 << "\n\n"; */

  // first: just use the scheme's in built function, being careful to use the correct logq
  e1 = scheme.ZZX_oo_difference(prodb, cprodb, decrypted.logq);
  // second: get the max coefficient of the difference
  diff = prodb - cprodb;
  e2 = polynomial_max_coeff(diff);
  // third: balance the difference and get the max coefficient
  diffb = scheme.balance_polynomial(diff, decrypted.logq);
  e3 = polynomial_max_coeff(diffb);
  // fourth: take the difference mod q instead
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq, prodb, cprodb, q, context.N);
  e4 = polynomial_max_coeff(diffq);
  /*cout << "MULT WITH BALANCING:\n";
  cout << "noise e1 = " << e1 << "\n";
  cout << "noise e2 = " << e2 << "\n";
  cout << "noise e3 = " << e3 << "\n";
  cout << "noise e4 = " << e4 << "\n\n"; */

  // here we start rescaling. Not needed for add noise!
  // now rescale both
  Ciphertext c_rs = scheme.reScaleBy(c3, logp);
  ZZX p_rs;
  Ring2Utils::rightShift(p_rs, prod, logp, context.N);
  // decrypt & print
  // cout << "Printing the rescale of r1*r2 mod Q:\n";
  // print_poly(p_rs, 12);
  // cout << "Printing the rescale of c1*c2:\n";
  decrypted = scheme.decrypt_no_decode(sk, c_rs);
  // print_poly(decrypted.mx, 12);

  // balance again
  prodb = scheme.balance_polynomial(p_rs, decrypted.logq);
  cprodb = scheme.balance_polynomial(decrypted.mx, decrypted.logq);

  // find noise now, using all 4 methods above
  // first: just use the scheme's in built function, being careful to use the correct logq
  e1 = scheme.ZZX_oo_difference(p_rs, decrypted.mx, decrypted.logq);
  // second: get the max coefficient of the difference
  diff = p_rs - decrypted.mx;
  e2 = polynomial_max_coeff(diff);
  // third: balance the difference and get the max coefficient
  diffb = scheme.balance_polynomial(diff, decrypted.logq);
  e3 = polynomial_max_coeff(diffb);
  // fourth: take the difference mod q instead
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq, p_rs, decrypted.mx, q, context.N);
  e4 = polynomial_max_coeff(diffq);

  /*cout << "\n######## NOISE AFTER RESCALE ########\n\n";
  cout << "WITHOUT BALANCING:\n";
  cout << "noise e1 = " << e1 << "\n";
  cout << "noise e2 = " << e2 << "\n";
  cout << "noise e3 = " << e3 << "\n";
  cout << "noise e4 = " << e4 << "\n\n"; */
  e1 = scheme.ZZX_oo_difference(prodb, cprodb, decrypted.logq);
  // second: get the max coefficient of the difference
  diff = prodb - cprodb;
  e2 = polynomial_max_coeff(diff);
  // third: balance the difference and get the max coefficient
  diffb = scheme.balance_polynomial(diff, decrypted.logq);
  e3 = polynomial_max_coeff(diffb);
  // fourth: take the difference mod q instead
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq, prodb, cprodb, q, context.N);
  e4 = polynomial_max_coeff(diffq);
  /*cout << "WITH BALANCING:\n";
  cout << "noise e1 = " << e1 << "\n";
  cout << "noise e2 = " << e2 << "\n";
  cout << "noise e3 = " << e3 << "\n";
  cout << "noise e4 = " << e4 << "\n\n"; */

  noises[4] = e1;
  noises[5] = e2;
  noises[6] = e3;
  noises[7] = e4;

  return noises;
}

vector<ZZ> noise_final_random_plaintexts(long logN, long logQ, long logp, ZZ plaintext_bound, SecretKey sk)
{
  // cout << "######## PARAMETERS: logN = " << logN << " logQ = " << logQ << " logDelta = " << logp << " bound = " << bound << "\n\n";
  Context context(logN, logQ);
  // SecretKey sk;
  Scheme scheme(sk, context);
  long slots = (1 << (logN - 1));

  // generate the random vecs that will encode to uniform random plaintexts
  complex<double> *r1 = generate_random_plaintext_v2(logN, plaintext_bound, logp, logQ);
  complex<double> *r2 = generate_random_plaintext_v2(logN, plaintext_bound, logp, logQ);
  complex<double> *r3 = generate_random_plaintext_v2(logN, plaintext_bound, logp, logQ);

  // Now encode for the plaintext circuit
  Plaintext m1 = context.encode(r1, slots, logp);
  Plaintext m2 = context.encode(r2, slots, logp);
  Plaintext m3 = context.encode(r3, slots, logp);
  // ZZX p_plain = m1.mx;
  // cout << "Printing p.mx\n\n";
  // print_poly(p_plain,12);

  Ciphertext c0 = scheme.encrypt(r1, slots, logp, logQ);
  Ciphertext c1 = scheme.encrypt(r2, slots, logp, logQ);
  Ciphertext c2 = scheme.encrypt(r3, slots, logp, logQ);
  // m1 = scheme.decrypt_no_decode(sk,c0);
  // p_plain = scheme.balance_polynomial(m1.mx,logQ);
  // cout << "Printing p.mx\n\n";
  // print_poly(p_plain,12);
  // Trying the vector
  vector<ZZ> noises(8);

  // Now add, and mult
  // logQ watch: keyswitch done modulo qQ, actual "multiplication" done modulo q (the "current" ciphertext modulus)

  Ciphertext c_add = scheme.add(c0, c1);
  Ciphertext c3 = scheme.mult(c2, c_add);
  // logQ watch: no reduction here

  // change: I'm going to do this modulo Q, the ciphertext modulus, to match the multiplication within scheme above
  ZZ q = context.qpowvec[c1.logq];
  ZZX add;  // add of the plaintexts
  ZZX prod; // prod of the plaintexts
  Ring2Utils::add(add, m1.mx, m2.mx, q, context.N);
  Ring2Utils::mult(prod, m3.mx, add, q, context.N);

  // printing add; no need for rescale, this shouldn't need anything else
  // cout << "printing r0 + r1 mod Q:\n\n";
  // print_poly(add, 12);
  Plaintext addp = add;

  // printing prod: this hasn't been rescaled yet
  // cout << "printing r1*r2 mod Q:\n\n";
  // print_poly(prod, 12);
  Plaintext prodp = prod;

  // running a decryption on c3, which is c_add*c2 without a rescale
  // also running a decryption of c_add = c_0 + c_1

  Plaintext decrypted_add = scheme.decrypt_no_decode(sk, c_add);
  Plaintext decrypted = scheme.decrypt_no_decode(sk, c3);

  // balance the addition (plaintext)
  ZZX addb = scheme.balance_polynomial(add, decrypted_add.logq);
  // balance ciphertext add balanced (ciphertext)
  ZZX caddb = scheme.balance_polynomial(decrypted_add.mx, decrypted_add.logq);

  // balance the prod (plaintext)
  ZZX prodb = scheme.balance_polynomial(prod, decrypted.logq);
  // balance ciphertext prod balanced (ciphertext)
  ZZX cprodb = scheme.balance_polynomial(decrypted.mx, decrypted.logq);

  /*
  cout << "Printing the decryption of c0+c1:\n\n";
  print_poly(addb, 12);
  cout << "Printing the plaintext addition of m1 + m2:\n\n";
  print_poly(caddb, 12);
  */

  // using fi for add noise
  // ei for mult noise

  // first: just use the scheme's in built function, being careful to use the correct logq
  ZZ f1 = scheme.ZZX_oo_difference(add, decrypted_add.mx, decrypted_add.logq);
  // second: get the max coefficient of the difference
  ZZX diff_add = add - decrypted_add.mx;
  ZZ f2 = polynomial_max_coeff(diff_add);
  // third: balance the difference and get the max coefficient
  ZZX diffb_add = scheme.balance_polynomial(diff_add, decrypted_add.logq);
  ZZ f3 = polynomial_max_coeff(diffb_add);
  // fourth: take the difference mod q instead
  ZZX diffq_add;
  q = context.qpowvec[decrypted_add.logq];
  Ring2Utils::sub(diffq_add, add, decrypted_add.mx, q, context.N);
  ZZ f4 = polynomial_max_coeff(diffq_add);
  /*cout << "\n######## ADD NOISE ########\n\n";
  cout << "WITHOUT BALANCING:\n";
  cout << "noise f1 = " << f1 << "\n";
  cout << "noise f2 = " << f2 << "\n";
  cout << "noise f3 = " << f3 << "\n";
  cout << "noise f4 = " << f4 << "\n\n"; */

  // first: just use the scheme's in built function, being careful to use the correct logq
  f1 = scheme.ZZX_oo_difference(addb, caddb, decrypted_add.logq);
  // second: get the max coefficient of the difference
  diff_add = addb - caddb;
  f2 = polynomial_max_coeff(diff_add);
  // third: balance the difference and get the max coefficient
  diffb_add = scheme.balance_polynomial(diff_add, decrypted_add.logq);
  f3 = polynomial_max_coeff(diffb_add);
  // fourth: take the difference mod q instead
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq_add, addb, caddb, q, context.N);
  f4 = polynomial_max_coeff(diffq_add);
  /*cout << "ADD WITH BALANCING:\n";
  cout << "noise f1 = " << f1 << "\n";
  cout << "noise f2 = " << f2 << "\n";
  cout << "noise f3 = " << f3 << "\n";
  cout << "noise f4 = " << f4 << "\n\n"; */

  noises[0] = f1;
  noises[1] = f2;
  noises[2] = f3;
  noises[3] = f4;

  // first: just use the scheme's in built function, being careful to use the correct logq
  ZZ e1 = scheme.ZZX_oo_difference(prod, decrypted.mx, decrypted.logq);
  // second: get the max coefficient of the difference
  ZZX diff = prod - decrypted.mx;
  ZZ e2 = polynomial_max_coeff(diff);
  // third: balance the difference and get the max coefficient
  ZZX diffb = scheme.balance_polynomial(diff, decrypted.logq);
  ZZ e3 = polynomial_max_coeff(diffb);
  // fourth: take the difference mod q instead
  ZZX diffq;
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq, prod, decrypted.mx, q, context.N);
  ZZ e4 = polynomial_max_coeff(diffq);
  /*cout << "\n######## MULT NOISE BEFORE RESCALE ########\n\n";
  cout << "WITHOUT BALANCING:\n";
  cout << "noise e1 = " << e1 << "\n";
  cout << "noise e2 = " << e2 << "\n";
  cout << "noise e3 = " << e3 << "\n";
  cout << "noise e4 = " << e4 << "\n\n"; */

  // first: just use the scheme's in built function, being careful to use the correct logq
  e1 = scheme.ZZX_oo_difference(prodb, cprodb, decrypted.logq);
  // second: get the max coefficient of the difference
  diff = prodb - cprodb;
  e2 = polynomial_max_coeff(diff);
  // third: balance the difference and get the max coefficient
  diffb = scheme.balance_polynomial(diff, decrypted.logq);
  e3 = polynomial_max_coeff(diffb);
  // fourth: take the difference mod q instead
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq, prodb, cprodb, q, context.N);
  e4 = polynomial_max_coeff(diffq);
  /*cout << "MULT WITH BALANCING:\n";
  cout << "noise e1 = " << e1 << "\n";
  cout << "noise e2 = " << e2 << "\n";
  cout << "noise e3 = " << e3 << "\n";
  cout << "noise e4 = " << e4 << "\n\n"; */

  // here we start rescaling. Not needed for add noise!
  // now rescale both
  Ciphertext c_rs = scheme.reScaleBy(c3, logp);
  ZZX p_rs;
  Ring2Utils::rightShift(p_rs, prod, logp, context.N);
  // decrypt & print
  // cout << "Printing the rescale of r1*r2 mod Q:\n";
  // print_poly(p_rs, 12);
  // cout << "Printing the rescale of c1*c2:\n";
  decrypted = scheme.decrypt_no_decode(sk, c_rs);
  // print_poly(decrypted.mx, 12);

  // balance again
  prodb = scheme.balance_polynomial(p_rs, decrypted.logq);
  cprodb = scheme.balance_polynomial(decrypted.mx, decrypted.logq);

  // find noise now, using all 4 methods above
  // first: just use the scheme's in built function, being careful to use the correct logq
  e1 = scheme.ZZX_oo_difference(p_rs, decrypted.mx, decrypted.logq);
  // second: get the max coefficient of the difference
  diff = p_rs - decrypted.mx;
  e2 = polynomial_max_coeff(diff);
  // third: balance the difference and get the max coefficient
  diffb = scheme.balance_polynomial(diff, decrypted.logq);
  e3 = polynomial_max_coeff(diffb);
  // fourth: take the difference mod q instead
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq, p_rs, decrypted.mx, q, context.N);
  e4 = polynomial_max_coeff(diffq);

  /*cout << "\n######## NOISE AFTER RESCALE ########\n\n";
  cout << "WITHOUT BALANCING:\n";
  cout << "noise e1 = " << e1 << "\n";
  cout << "noise e2 = " << e2 << "\n";
  cout << "noise e3 = " << e3 << "\n";
  cout << "noise e4 = " << e4 << "\n\n"; */
  e1 = scheme.ZZX_oo_difference(prodb, cprodb, decrypted.logq);
  // second: get the max coefficient of the difference
  diff = prodb - cprodb;
  e2 = polynomial_max_coeff(diff);
  // third: balance the difference and get the max coefficient
  diffb = scheme.balance_polynomial(diff, decrypted.logq);
  e3 = polynomial_max_coeff(diffb);
  // fourth: take the difference mod q instead
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq, prodb, cprodb, q, context.N);
  e4 = polynomial_max_coeff(diffq);
  /*cout << "WITH BALANCING:\n";
  cout << "noise e1 = " << e1 << "\n";
  cout << "noise e2 = " << e2 << "\n";
  cout << "noise e3 = " << e3 << "\n";
  cout << "noise e4 = " << e4 << "\n\n"; */

  noises[4] = e1;
  noises[5] = e2;
  noises[6] = e3;
  noises[7] = e4;

  return noises;
}

vector<ZZ> noise_final_fixed_vector(long logN, long logQ, long logp, SecretKey sk)
{
  // cout << "######## PARAMETERS: logN = " << logN << " logQ = " << logQ << " logDelta = " << logp << " bound = " << bound << "\n\n";
  Context context(logN, logQ);
  // SecretKey sk;
  Scheme scheme(sk, context);
  long slots = (1 << (logN - 1));

  // Now encode
  // Plaintext r1 = attained_polynomial(logN, bound, logp, logQ);
  // Plaintext r2 = attained_polynomial(logN, bound, logp, logQ);
  // Plaintext r3 = attained_polynomial(logN, bound, logp, logQ);
  complex<double> *vec = constant_poly_vector(logN);

  Plaintext r1 = context.encode(vec, slots, logp);
  Plaintext r2 = context.encode(vec, slots, logp);
  Plaintext r3 = context.encode(vec, slots, logp);
  // print_poly(r1.mx,10);

  // Now encrypt
  // I know the subscripts are a bit weird. But I wanted
  // to keep c3 as in new-stuff.h to avoid potential errors
  Ciphertext c0 = scheme.encrypt(vec, slots, logp, logQ);
  Ciphertext c1 = scheme.encrypt(vec, slots, logp, logQ);
  Ciphertext c2 = scheme.encrypt(vec, slots, logp, logQ);

  // Trying the vector
  vector<ZZ> noises(8);

  // Now add, and mult
  // logQ watch: keyswitch done modulo qQ, actual "multiplication" done modulo q (the "current" ciphertext modulus)

  Ciphertext c_add = scheme.add(c0, c1);
  Ciphertext c3 = scheme.mult(c2, c_add);
  // logQ watch: no reduction here

  // change: I'm going to do this modulo Q, the ciphertext modulus, to match the multiplication within scheme above
  ZZ q = context.qpowvec[c1.logq];
  ZZX add;  // add of the plaintexts
  ZZX prod; // prod of the plaintexts
  Ring2Utils::add(add, r1.mx, r2.mx, q, context.N);
  Ring2Utils::mult(prod, r3.mx, add, q, context.N);

  // cout << "printing r0 + r1 mod Q:\n\n";
  // print_poly(add, 12);
  Plaintext addp = add;

  // printing prod: this hasn't been rescaled yet
  // cout << "printing r1*r2 mod Q:\n\n";
  // print_poly(prod, 12);
  Plaintext prodp = prod;

  // running a decryption on c3, which is c_add*c2 without a rescale
  // also running a decryption of c_add = c_0 + c_1

  Plaintext decrypted_add = scheme.decrypt_no_decode(sk, c_add);
  Plaintext decrypted = scheme.decrypt_no_decode(sk, c3);

  // cout << "Printing the decryption of c0+c1:\n\n";
  // print_poly(decrypted_add.mx, 12);
  // cout << "Printing the decryption of c_add*c2:\n\n";
  // print_poly(decrypted.mx, 12);

  // balance the addition (plaintext)
  ZZX addb = scheme.balance_polynomial(add, decrypted_add.logq);
  // balance ciphertext add balanced (ciphertext)
  ZZX caddb = scheme.balance_polynomial(decrypted_add.mx, decrypted_add.logq);

  // balance the prod (plaintext)
  ZZX prodb = scheme.balance_polynomial(prod, decrypted.logq);
  // balance ciphertext prod balanced (ciphertext)
  ZZX cprodb = scheme.balance_polynomial(decrypted.mx, decrypted.logq);

  // now manually get the noise. Do this a variety of ways and see what shakes out
  // I think all of these are the same for small enough noise but potentially will be different once
  // the noise wraps around q/2....

  // using fi for add noise
  // ei for mult noise

  // first: just use the scheme's in built function, being careful to use the correct logq
  ZZ f1 = scheme.ZZX_oo_difference(add, decrypted_add.mx, decrypted_add.logq);
  // second: get the max coefficient of the difference
  ZZX diff_add = add - decrypted_add.mx;
  ZZ f2 = polynomial_max_coeff(diff_add);
  // third: balance the difference and get the max coefficient
  ZZX diffb_add = scheme.balance_polynomial(diff_add, decrypted_add.logq);
  ZZ f3 = polynomial_max_coeff(diffb_add);
  // fourth: take the difference mod q instead
  ZZX diffq_add;
  q = context.qpowvec[decrypted_add.logq];
  Ring2Utils::sub(diffq_add, add, decrypted_add.mx, q, context.N);
  ZZ f4 = polynomial_max_coeff(diffq_add);
  /*cout << "\n######## ADD NOISE ########\n\n";
  cout << "WITHOUT BALANCING:\n";
  cout << "noise f1 = " << f1 << "\n";
  cout << "noise f2 = " << f2 << "\n";
  cout << "noise f3 = " << f3 << "\n";
  cout << "noise f4 = " << f4 << "\n\n"; */

  // first: just use the scheme's in built function, being careful to use the correct logq
  f1 = scheme.ZZX_oo_difference(addb, caddb, decrypted_add.logq);
  // second: get the max coefficient of the difference
  diff_add = addb - caddb;
  f2 = polynomial_max_coeff(diff_add);
  // third: balance the difference and get the max coefficient
  diffb_add = scheme.balance_polynomial(diff_add, decrypted_add.logq);
  f3 = polynomial_max_coeff(diffb_add);
  // fourth: take the difference mod q instead
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq_add, addb, caddb, q, context.N);
  f4 = polynomial_max_coeff(diffq_add);
  /*cout << "ADD WITH BALANCING:\n";
  cout << "noise f1 = " << f1 << "\n";
  cout << "noise f2 = " << f2 << "\n";
  cout << "noise f3 = " << f3 << "\n";
  cout << "noise f4 = " << f4 << "\n\n"; */

  noises[0] = f1;
  noises[1] = f2;
  noises[2] = f3;
  noises[3] = f4;

  // first: just use the scheme's in built function, being careful to use the correct logq
  ZZ e1 = scheme.ZZX_oo_difference(prod, decrypted.mx, decrypted.logq);
  // second: get the max coefficient of the difference
  ZZX diff = prod - decrypted.mx;
  ZZ e2 = polynomial_max_coeff(diff);
  // third: balance the difference and get the max coefficient
  ZZX diffb = scheme.balance_polynomial(diff, decrypted.logq);
  ZZ e3 = polynomial_max_coeff(diffb);
  // fourth: take the difference mod q instead
  ZZX diffq;
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq, prod, decrypted.mx, q, context.N);
  ZZ e4 = polynomial_max_coeff(diffq);
  /*cout << "\n######## MULT NOISE BEFORE RESCALE ########\n\n";
  cout << "WITHOUT BALANCING:\n";
  cout << "noise e1 = " << e1 << "\n";
  cout << "noise e2 = " << e2 << "\n";
  cout << "noise e3 = " << e3 << "\n";
  cout << "noise e4 = " << e4 << "\n\n"; */

  // first: just use the scheme's in built function, being careful to use the correct logq
  e1 = scheme.ZZX_oo_difference(prodb, cprodb, decrypted.logq);
  // second: get the max coefficient of the difference
  diff = prodb - cprodb;
  e2 = polynomial_max_coeff(diff);
  // third: balance the difference and get the max coefficient
  diffb = scheme.balance_polynomial(diff, decrypted.logq);
  e3 = polynomial_max_coeff(diffb);
  // fourth: take the difference mod q instead
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq, prodb, cprodb, q, context.N);
  e4 = polynomial_max_coeff(diffq);
  /*cout << "MULT WITH BALANCING:\n";
  cout << "noise e1 = " << e1 << "\n";
  cout << "noise e2 = " << e2 << "\n";
  cout << "noise e3 = " << e3 << "\n";
  cout << "noise e4 = " << e4 << "\n\n"; */

  // here we start rescaling. Not needed for add noise!
  // now rescale both
  Ciphertext c_rs = scheme.reScaleBy(c3, logp);
  ZZX p_rs;
  Ring2Utils::rightShift(p_rs, prod, logp, context.N);
  // decrypt & print
  // cout << "Printing the rescale of r1*r2 mod Q:\n";
  // print_poly(p_rs, 12);
  // cout << "Printing the rescale of c1*c2:\n";
  decrypted = scheme.decrypt_no_decode(sk, c_rs);
  // print_poly(decrypted.mx, 12);

  // balance again
  prodb = scheme.balance_polynomial(p_rs, decrypted.logq);
  cprodb = scheme.balance_polynomial(decrypted.mx, decrypted.logq);

  // find noise now, using all 4 methods above
  // first: just use the scheme's in built function, being careful to use the correct logq
  e1 = scheme.ZZX_oo_difference(p_rs, decrypted.mx, decrypted.logq);
  // second: get the max coefficient of the difference
  diff = p_rs - decrypted.mx;
  e2 = polynomial_max_coeff(diff);
  // third: balance the difference and get the max coefficient
  diffb = scheme.balance_polynomial(diff, decrypted.logq);
  e3 = polynomial_max_coeff(diffb);
  // fourth: take the difference mod q instead
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq, p_rs, decrypted.mx, q, context.N);
  e4 = polynomial_max_coeff(diffq);

  /*cout << "\n######## NOISE AFTER RESCALE ########\n\n";
  cout << "WITHOUT BALANCING:\n";
  cout << "noise e1 = " << e1 << "\n";
  cout << "noise e2 = " << e2 << "\n";
  cout << "noise e3 = " << e3 << "\n";
  cout << "noise e4 = " << e4 << "\n\n"; */
  e1 = scheme.ZZX_oo_difference(prodb, cprodb, decrypted.logq);
  // second: get the max coefficient of the difference
  diff = prodb - cprodb;
  e2 = polynomial_max_coeff(diff);
  // third: balance the difference and get the max coefficient
  diffb = scheme.balance_polynomial(diff, decrypted.logq);
  e3 = polynomial_max_coeff(diffb);
  // fourth: take the difference mod q instead
  q = context.qpowvec[decrypted.logq];
  Ring2Utils::sub(diffq, prodb, cprodb, q, context.N);
  e4 = polynomial_max_coeff(diffq);
  /*cout << "WITH BALANCING:\n";
  cout << "noise e1 = " << e1 << "\n";
  cout << "noise e2 = " << e2 << "\n";
  cout << "noise e3 = " << e3 << "\n";
  cout << "noise e4 = " << e4 << "\n\n"; */

  noises[4] = e1;
  noises[5] = e2;
  noises[6] = e3;
  noises[7] = e4;

  delete[] vec;
  return noises;
}

vector<double> complex_to_complex(long logN, long logQ, long logp, double bound = 1.0)
{
  Context context(logN, logQ);
  SecretKey sk;
  Scheme scheme(sk, context);
  long slots = (1 << (logN - 1));

  // Generate some random (real) data to be operated on
  complex<double> *m1 = EvaluatorUtils::randomComplexArray(slots, bound);
  complex<double> *m2 = EvaluatorUtils::randomComplexArray(slots, bound);
  complex<double> *m3 = EvaluatorUtils::randomComplexArray(slots, bound);

  // Now encrypt
  // I know the subscripts are a bit weird. But I wanted
  // to keep c3 as in new-stuff.h to avoid potential errors
  Ciphertext c1 = scheme.encrypt(m1, slots, logp, logQ);
  Ciphertext c2 = scheme.encrypt(m2, slots, logp, logQ);
  Ciphertext c3 = scheme.encrypt(m3, slots, logp, logQ);

  // Trying the vector
  vector<double> precision(4);

  // Now let's operate on raw data
  complex<double> *add_raw = complex_vector_add(m1, m2, slots);
  complex<double> *mult_raw = complex_vector_mult(m3, add_raw, slots);

  // Now add, and mult
  // logQ watch: keyswitch done modulo qQ, actual "multiplication" done modulo q (the "current" ciphertext modulus)

  Ciphertext c_add = scheme.add(c1, c2);
  Ciphertext c_mult_no_rescale = scheme.mult(c3, c_add);

  // let's rescale now
  Ciphertext c_mult = scheme.reScaleBy(c_mult_no_rescale, logp);

  // can't remember what this did, let's try to comment it out
  // ZZ q = context.qpowvec[c1.logq];

  // decrypting
  complex<double> *decrypted_add = scheme.decrypt(sk, c_add);
  complex<double> *decrypted_mult = scheme.decrypt(sk, c_mult);

  // now let's get the precision loss, hopefully
  double precision_loss_difference_add = max_difference(add_raw, decrypted_add, slots);
  double precision_loss_real_difference_add = max_real_difference(add_raw, decrypted_add, slots);

  double precision_loss_difference_mult = max_difference(mult_raw, decrypted_mult, slots);
  double precision_loss_real_difference_mult = max_real_difference(add_raw, decrypted_mult, slots);

  precision[0] = precision_loss_difference_add;
  precision[1] = precision_loss_real_difference_add;
  precision[2] = precision_loss_difference_mult;
  precision[3] = precision_loss_real_difference_add;

  return precision;
}

vector<ZZ> noise_fresh_random_plaintexts(long logN, long logQ, long logp, ZZ plaintext_bound, SecretKey sk)
{
  // cout << "######## PARAMETERS: logN = " << logN << " logQ = " << logQ << " logDelta = " << logp << " bound = " << bound << "\n\n";
  Context context(logN, logQ);
  // SecretKey sk;
  Scheme scheme(sk, context);
  long slots = (1 << (logN - 1));

  // Generate random plaintexts
  // Plaintext m1 = generate_random_plaintext(logN, plaintext_bound, logp, logQ);
  complex<double> *vec = generate_random_plaintext_v2(logN, plaintext_bound, logp, logQ);
  Plaintext p = context.encode(vec, slots, logp);
  ZZX p_plain = p.mx;
  cout << "Printing p.mx\n\n";
  print_poly(p_plain, 12);
  // Now encrypt
  // I know the subscripts are a bit weird. But I wanted
  // to keep c3 as in new-stuff.h to avoid potential errors
  Ciphertext c1 = scheme.encrypt(vec, slots, logp, logQ);
  Ciphertext c2 = c1;
  ZZ q = context.qpowvec[c1.logq];
  ZZ qq = context.qpowvec[c2.logq];

  // Trying the vector
  vector<ZZ> noises(8);

  // running a decryption on c3, which is c_add*c2 without a rescale
  // also running a decryption of c_add = c_0 + c_1

  Plaintext decrypted_fresh = scheme.decrypt_no_decode(sk, c1);
  Plaintext decrypted_test = scheme.decryptMsg(sk, c2);
  // cout << "printing decrypted_fresh.mx:\n\n";
  // print_poly(decrypted_fresh.mx,10);
  // cout << "\nprinting decrypted_test.mx:\n\n";
  // print_poly(decrypted_test.mx,10);
  ZZX freshb = scheme.balance_polynomial(decrypted_fresh.mx, c1.logq);
  cout << "\nprinting decrypted_fresh.mx after balancing:\n\n";
  print_poly(freshb, 10);
  // now manually get the noise. Do this a variety of ways and see what shakes out
  // I think all of these are the same for small enough noise but potentially will be different once
  // the noise wraps around q/2....

  // using fi for add noise
  // ei for mult noise

  // first: just use the scheme's in built function, being careful to use the correct logq
  ZZ f1 = scheme.ZZX_oo_difference(p_plain, freshb, decrypted_fresh.logq);
  // second: get the max coefficient of the difference
  ZZX diff = p_plain - freshb;
  ZZ f2 = polynomial_max_coeff(diff);
  // third: balance the difference and get the max coefficient
  ZZX diffb = scheme.balance_polynomial(diff, decrypted_fresh.logq);
  ZZ f3 = polynomial_max_coeff(diffb);
  // fourth: take the difference mod q instead
  ZZX diffq;
  q = context.qpowvec[decrypted_fresh.logq];
  Ring2Utils::sub(diffq, p_plain, freshb, q, context.N);
  ZZ f4 = polynomial_max_coeff(diffq);
  cout << "\n######## FRESH NOISE ########\n\n";
  cout << "WITHOUT BALANCING:\n";
  cout << "noise f1 = " << f1 << "\n";
  cout << "noise f2 = " << f2 << "\n";
  cout << "noise f3 = " << f3 << "\n";
  cout << "noise f4 = " << f4 << "\n\n";

  noises[0] = f1;
  noises[1] = f2;
  noises[2] = f3;
  noises[3] = f4;

  // first: just use the scheme's in built function, being careful to use the correct logq
  ZZ e1 = scheme.ZZX_oo_difference(p_plain, freshb, decrypted_test.logq);
  // second: get the max coefficient of the difference
  ZZX diff_test = p_plain - freshb;
  ZZ e2 = polynomial_max_coeff(diff_test);
  // third: balance the difference and get the max coefficient
  ZZX diffb_test = scheme.balance_polynomial(diff_test, decrypted_test.logq);
  ZZ e3 = polynomial_max_coeff(diffb_test);
  // fourth: take the difference mod q instead
  ZZX diffq_test;
  Ring2Utils::sub(diffq_test, p_plain, freshb, qq, context.N);
  ZZ e4 = polynomial_max_coeff(diffq_test);

  noises[4] = e1;
  noises[5] = e2;
  noises[6] = e3;
  noises[7] = e4;
  cout << "\n######## FRESH NOISE ########\n\n";
  cout << "WITH BALANCING:\n";
  cout << "noise e1 = " << f1 << "\n";
  cout << "noise e2 = " << f2 << "\n";
  cout << "noise e3 = " << f3 << "\n";
  cout << "noise e4 = " << f4 << "\n\n";

  delete[] vec;
  return noises;
}

void noise_final_average(long logN, long logQ, long logp, ZZ plaintext_bound, int loop, SecretKey sk, int fn)
{
  // now let us loop
  // only doing after rescale but this can be changed

  // the int fn should be one of 0,1,2, where 0 calls fixed_vector, 1 calls random_plaintext, 2 calls random_reals

  vector<ZZ> results(8);
  vector<ZZ> max_values(8);
  for (int i = 0; i < loop; i++)
  {
    vector<ZZ> noise;
    switch (fn)
    {
    case 0:
      noise = noise_final_fixed_vector(logN, logQ, logp, sk);
      break;
    case 1:
      noise = noise_final_random_plaintexts(logN, logQ, logp, plaintext_bound, sk);
      break;
      // case 2:
      // noise = noise_final_random_reals(logN, logQ, logp, bound,sk);
      // break;
      // case 3:
      // noise = noise_fresh_random_plaintexts(logN, logQ, logp, plaintext_bound,sk);
      // default:
      return;
    }
    results[0] += noise[0];
    results[1] += noise[1];
    results[2] += noise[2];
    results[3] += noise[3];
    results[4] += noise[4];
    results[5] += noise[5];
    results[6] += noise[6];
    results[7] += noise[7];

    for (int j = 0; j < 8; j++)
    {
      if (noise[j] > max_values[j])
      {
        max_values[j] = noise[j];
      }
    }
  }
  for (int j = 0; j < 8; j++)
  {
    results[j] = results[j] / loop;
  }

  cout << "The average noise values for addition are:\n";
  cout << "noise f1 = " << results[0] << "\n";
  cout << "noise f2 = " << results[1] << "\n";
  cout << "noise f3 = " << results[2] << "\n";
  cout << "noise f4 = " << results[3] << "\n\n";

  cout << "The max noise values for addition are:\n";
  cout << "noise f1 = " << max_values[0] << "\n";
  cout << "noise f2 = " << max_values[1] << "\n";
  cout << "noise f3 = " << max_values[2] << "\n";
  cout << "noise f4 = " << max_values[3] << "\n\n";

  cout << "The average noise values for multiplication are:\n";
  cout << "noise e1 = " << results[4] << "\n";
  cout << "noise e2 = " << results[5] << "\n";
  cout << "noise e3 = " << results[6] << "\n";
  cout << "noise e4 = " << results[7] << "\n\n";

  cout << "The max noise values for multiplication are:\n";
  cout << "noise e1 = " << max_values[4] << "\n";
  cout << "noise e2 = " << max_values[5] << "\n";
  cout << "noise e3 = " << max_values[6] << "\n";
  cout << "noise e4 = " << max_values[7] << "\n\n";
}
