# This file contains the functions to derive the prior canonincal embedding bounds of Table 6 (eprint version) for Textbook CKKS (HEAAN) and Table 7 (eprint version) for FullRNS-HEAAN
# These bounds from prior work can be compared to our bounds
from numpy import *
from math import *

# The fresh bound is the same for Textbook CKKS (HEAAN) and FullRNS-HEAAN
def fresh_heaan(sigma, N, h):
    return 8 * sqrt(2) * sigma * N + 6 * sigma * sqrt(N) + 16 * sigma * sqrt(h * N)

# The add bound is the same for Textbook CKKS (HEAAN) and FullRNS-HEAAN
def add(input1, input2):
    return input1 + input2

# The pre_mult bound is the same for Textbook CKKS (HEAAN) and FullRNS-HEAAN
def pre_mult(input1, input2, nu1, nu2):
    return nu1 * input2 + nu2 * input1 + input1 * input2

# Textbook CKKS (HEAAN) key switch
def ks_textbook(input1, P, N, h, sigma, Q_ell, round=False):
    ks_additive_noise_bound = (8 / sqrt(3)) * N * sigma * Q_ell/P
    if round == True:
        ks_additive_noise_bound += sqrt(3 * N) + (8 / sqrt(3)) * sqrt(h * N)
    return input1 + ks_additive_noise_bound

# FullRNS-HEAAN key switch
def ks_fullrns_heaan(input1, k, N, h, L, Q_L, P, sigma):
    return input1 + (8 / sqrt(3)) * N * sigma * Q_L/P * sqrt(L * L + 1) + (k + 1) * (sqrt(3 * N) + (8 / sqrt(3)) * sqrt(h * N))

# Textbook CKKS (HEAAN) rescale
def rs_textbook(input1, delta, N, h):
    return input1/delta + sqrt(3 * N) + (8 / sqrt(3)) * sqrt(h * N)

# FullRNS-HEAAN rescale
# NB Q_ell = prod_{i=1 to q_ell} q_i for q_i the RNS moduli
def rs_fullrns_heaan(input1, q_ell, N, h):
    return input1/q_ell + sqrt(3 * N) + (8 / sqrt(3)) * sqrt(h * N)


# Top level function for Textbook HEAAN
def priorCE_textbook(delta, input_bound, N, sigma, h, Q_ell, tag=None):
    # For mult need bounds nu1 and nu2 on ||m||^can where m in the ring is the input message to mult. At fresh this should be delta * input_bound, after add this should be 2 * delta * input_bound
    fresh_message_bound = input_bound * delta

    # NB: prior work did not take into account encoding noise
    bound_fresh = fresh_heaan(N, sigma, h)
    bound_fresh = rs_textbook(bound_fresh, Q_ell, N, h)
    bound_add = add(bound_fresh, bound_fresh)
    bound_premult = pre_mult(bound_fresh, bound_add, fresh_message_bound, sqrt(
        2) * fresh_message_bound)  # ||a + b || < sqrt(||a||^2 + ||b||^2)
    # P = Q_L for HEAAN and we are calling ks at ell = L
    bound_ks = ks_textbook(bound_premult, Q_ell, N, h, sigma, Q_ell)
    bound_rs = rs_textbook(bound_ks, delta, N, h)

    decode_add = log(bound_add, 2) - log(delta, 2)
    decode_mult = log(bound_rs, 2) - log(delta, 2)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N, 2)
    result["log(Q_ell)"] = round(log(Q_ell, 2), 1)
    result["log(Δ)"] = log(delta, 2)
    result["log(MB)"] = log(fresh_message_bound, 2)
    result["r+"] = decode_add
    result["rx"] = decode_mult
    result["c+"] = decode_add
    result["cx"] = decode_mult
    result["Rq+"] = round(log(bound_add, 2), 2)
    result["Rqx"] = round(log(bound_rs, 2), 2)

    return result


# Top level function for FullRNS-HEAAN
def priorCE_fullrns_heaan(delta, input_bound, N, sigma, h, q_ell, k, L, Q_L, P, tag=None):

    # For mult need bounds nu1 and nu2 on ||m||^can where m in the ring is the input message to mult. At fresh this should be delta * input_bound, after add this should be 2 * delta * input_bound
    fresh_message_bound = input_bound * delta

    # NB: prior work did not take into account encoding noise
    bound_fresh = fresh_heaan(N, sigma, h)
    bound_add = add(bound_fresh, bound_fresh)
    bound_premult = pre_mult(bound_fresh, bound_add,
                             fresh_message_bound, 2 * fresh_message_bound)
    bound_ks = ks_fullrns_heaan(bound_premult, k, N, h, L, Q_L, P, sigma)
    bound_rs = rs_fullrns_heaan(bound_ks, q_ell, N, h)

    decode_add = log(bound_add, 2) - log(delta, 2)
    decode_mult = log(bound_rs, 2) - log(delta, 2)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N, 2)
    result["log(q_ell)"] = round(log(q_ell, 2), 1)
    result["log(Δ)"] = log(delta, 2)
    result["log(MB)"] = log(fresh_message_bound, 2)
    result["c+"] = decode_add
    result["cx"] = decode_mult

    return result