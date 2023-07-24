# This file contains the functions to derive the CE bounds of Table 1 (eprint version) for Textbook CKKS (HEAAN), Table 3 (eprint version) for FullRNS-HEAAN

from math import *

# H_C(alpha, N) as defined above Theorem 2 (eprint version) 
def H_C(alpha, N):
    return sqrt(-log(1 - pow(1 - alpha, 2. / N)))

# Fresh bound is the same for Textbook (HEAAN) and FullRNS-HEAAN 
def fresh_heaan(sigma, N, alpha, h):
    return sigma * sqrt((1.0/2)*N**2 + h * N + N) * H_C(alpha, N)

# Add bound is the same for Textbook CKKS (HEAAN) and FullRNS-HEAAN
def add(input1, input2):
    return input1 + input2

# Premult bound is the same for Textbook CKKS (HEAAN) and FullRNS-HEAAN
# message_bound1 and message_bound2 are bounds on ||m1||^{can}_infty and ||m2||^{can}_infty
def premult(input1, input2, message_bound1, message_bound2):
    return input1 * message_bound2 + input2 * message_bound1 + input1 * input2

# In Textbook CKKS key switching as implemented in HEAAN we will not need rounding as P = Q_ell
def ks_textbook(input1, N, sigma, alpha, Q_ell, P, h, round=False):
    eta_ks_sq = (Q_ell / P) * (Q_ell / P) * N * sigma * sigma
    if round == True:
        eta_ks_sq += h + 1
    eta_ks_sq *= (1.0/12)
    eta_ks = sqrt(eta_ks_sq)
    return input1 + sqrt(N) * eta_ks * H_C(alpha, N)

# Key switching for RNS variant following Lemma 9 (eprint version) in the situation of FullRNS-HEAAN
def ks_fullrns_heaan(input1, k, N, h, alpha, L, Q_L, P, sigma): 
    eta_ks_sq = (Q_L / P) * (Q_L / P) * N * sigma * sigma * (L * L + 1)
    eta_ks_sq += (k * k + 1) * (h + 1)
    eta_ks_sq *= (1.0/12)
    eta_ks = sqrt(eta_ks_sq)
    return input1 + sqrt(N) * eta_ks * H_C(alpha, N)

# Textbook CKKS rescale as implemented in HEAAN 
def rs_textbook(input1, delta, N, h, alpha):
    return input1 / delta + sqrt(N * (1.0/12) * (h + 1)) * H_C(alpha, N)

# Key switching for RNS variant following Lemma 10 (eprint version) in the situation of FullRNS-HEAAN
def rs_fullrns_heaan(input1, q_ell, N, h, alpha):
    return input1 / q_ell + sqrt(N * (1.0/12) * (h + 1)) * H_C(alpha, N)


# Functions for determining encoding error
def can_emb_encode(N):
    return N / pi

def can_emb_decode(noise_bound, delta):
    return noise_bound / delta

# Top level function for Textbook HEAAN
def CE_textbook(N, B, sigma, delta, alpha, h, Q_ell, tag=None):
    # B is a bound on the input message in the canonical embedding, so bound ||m||^can is given by Delta * B
    fresh_message_bound = B * delta
    fresh_bound = fresh_heaan(sigma, N, alpha, h)
    fresh_bound += can_emb_encode(N)  # Account for encoding noise
    fresh_bound = rs_textbook(fresh_bound, Q_ell, N, h, alpha)
    add_bound = add(fresh_bound, fresh_bound)
    premult_bound = premult(fresh_bound, add_bound,
                            fresh_message_bound, 2 * fresh_message_bound) # ||a + b || < ||a|| + ||b||
    ks_bound = ks_textbook(premult_bound, N, sigma, alpha,
                           Q_ell, Q_ell, h, round=False)
    rs_bound = rs_textbook(ks_bound, delta, N, h, alpha)

    decode_add = can_emb_decode(add_bound, delta)
    decode_mult = can_emb_decode(rs_bound, delta)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N, 2)
    result["log(Δ)"] = log(delta, 2)
    result["α"] = alpha
    result["log(MB)"] = log(fresh_message_bound, 2)
    result["r+"] = round(log(decode_add, 2), 2)
    result["rx"] = round(log(decode_mult, 2), 2)
    result["c+"] = round(log(decode_add, 2), 2)
    result["cx"] = round(log(decode_mult, 2), 2)
    result["Rq+"]= round(log(add_bound, 2), 2)
    result["Rqx"]= round(log(rs_bound, 2), 2)

    return result

# Top level function for FullRNS-HEAAN
def CE_fullrns_heaan(N, B, sigma, delta, alpha, h, k, q_ell, L, Q_L, P, tag=None):    
    fresh_message_bound = B * delta
    fresh_bound = fresh_heaan(sigma, N, alpha, h)
    fresh_bound += can_emb_encode(N)  # Account for encoding noise
    add_bound = add(fresh_bound, fresh_bound)
    premult_bound = premult(fresh_bound, add_bound,
                            fresh_message_bound, 2 * fresh_message_bound)
    ks_bound = ks_fullrns_heaan(premult_bound, k, N, h, alpha, L, Q_L, P, sigma)
    rs_bound = rs_fullrns_heaan(ks_bound, q_ell, N, h, alpha)

    decode_add = can_emb_decode(add_bound, delta)
    decode_mult = can_emb_decode(rs_bound, delta)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N, 2)
    result["log(Δ)"] = log(delta, 2)
    result["α"] = alpha
    result["log(MB)"] = log(fresh_message_bound, 2)
    result["+"] = round(log(decode_add, 2), 2)
    result["x"] = round(log(decode_mult, 2), 2)

    return result