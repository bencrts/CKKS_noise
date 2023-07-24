# This file contains the functions to derive the WCR bounds of Table 1 (eprint version) for Textbook CKKS (HEAAN)
from math import *
from scipy.special import erfinv

# H_R(alpha, N) as defined above Theorem 2 (eprint version) 
def H_R(alpha, N):
    return erfinv(pow(1 - alpha, (1.0 / N)))

# Fresh bound
def fresh(sigma, N, h, alpha):
    return sigma * sqrt(N + 2 * h + 2) * H_R(alpha, N)
# Add bound
def add(input1, input2):
    return input1 + input2

# Premult bound
# message_bound1 and message_bound2 are bounds on ||m1||_infty and ||m2||_infty in the ring
def premult(input1, input2, message_bound1, message_bound2, N):
    return N * input1 * message_bound2 + N * input2 * message_bound1 + N * input1 * input2

# Key switch bound
def ks(input1, N, sigma, alpha, Q_ell, P, h, round=False):  # We will not need rounding as P = Q_ell
    eta_ks_sq = (Q_ell / P) * (Q_ell / P) * N * sigma * sigma
    if round == True:
        eta_ks_sq += h + 1
    eta_ks_sq *= (1.0/12)
    eta_ks = sqrt(eta_ks_sq)
    return input1 + sqrt(2) * eta_ks * H_R(alpha, N)

# Rescale bound
def rs(input1, delta, h, N, alpha):
    return input1 / delta + sqrt((1.0/6) * (h + 1)) * H_R(alpha, N)

# Will interpret the results in CC hence need to decode
def wcr_decode(N, noise_bound, delta):
    return (2 * N * noise_bound) / (pi * delta)

# Top level function for Textbook HEAAN
def WCR_textbook(N, B, sigma, delta, alpha, h, Q_ell, tag=None):
    # B is a bound in CC, and for premult we require bound on ring coeffs. We have ||m|| \leq ||m||_can \leq B * delta
    fresh_message_bound = B * delta
    fresh_bound = fresh(sigma, N, h, alpha)
    fresh_bound += 0.5  # Account for encoding noise
    fresh_bound = rs(fresh_bound, Q_ell, h, N, alpha)
    add_bound = add(fresh_bound, fresh_bound)
    premult_bound = premult(fresh_bound, add_bound,
                            fresh_message_bound, 2 * fresh_message_bound, N)
    ks_bound = ks(premult_bound, N, sigma, alpha, Q_ell, Q_ell, h, round=False)
    rs_bound = rs(ks_bound, delta, h, N, alpha)

    decode_add = wcr_decode(N, add_bound, delta)
    decode_rs = wcr_decode(N, rs_bound, delta)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N, 2)
    result["log(Δ)"] = log(delta, 2)
    result["α"] = alpha
    result["log(MB)"] = round(log(fresh_message_bound, 2), 2)
    result["r+"] = round(log(decode_add, 2), 2)
    result["rx"] = round(log(decode_rs, 2), 2)
    result["c+"] = round(log(decode_add, 2), 2)
    result["cx"] = round(log(decode_rs, 2), 2)
    result["Rq+"] = round(log(add_bound, 2), 2)
    result["Rqx"] = round(log(rs_bound, 2), 2)

    return result
