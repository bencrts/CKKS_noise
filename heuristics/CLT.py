from math import sqrt, log
from scipy.special import erfinv
from config import *

# This file contains the functions to derive the CLT noise analysis of Table 2 (eprint version) for Textbook CKKS (HEAAN) and Table 3 (eprint version) for FullRNS-HEAAN
# Each function returns a value for the variance rho^2.
# Once a final rho is calculated, this can be translated to a worst case bound on error coeffs, a total error in \CC^N/2, or a real error in \CC^N/2.

# The fresh variance is the same in Textbook CKKS (HEAAN) and FullRNS-HEAAN
def fresh_heaan(sigma, N, h):
    return sigma * sigma * ((1.0/2) * N + h + 1)

# The add variance is the same in Textbook CKKS (HEAAN) and FullRNS-HEAAN
def add(sq_rho_1, sq_rho_2):  # Requires independence
    return sq_rho_1 + sq_rho_2

# The pre_mult variance is the same in Textbook CKKS (HEAAN) and FullRNS-HEAAN
# M_1 and M_2 are (bounds on) |m_1|^2, |m_2|^2 the squared 2 norm of the polynomial
def mult(sq_rho_1, sq_rho_2, M_1, M_2, N):
    return N * sq_rho_1 * sq_rho_2 + sq_rho_1 * M_2 + sq_rho_2 * M_1

# In Textbook CKKS key switching as implemented in HEAAN we will not need rounding as P = Q_ell
def ks_textbook(sq_rho, N, h, sigma, Q_ell, P, round=False):
    eta_ks_sq = (Q_ell / P) * (Q_ell / P) * N * sigma * sigma
    if round == True:
        eta_ks_sq += h + 1
    eta_ks_sq *= (1.0/12)
    return sq_rho + eta_ks_sq

# Key switching for RNS variant following Lemma 9 (eprint version) in the situation of FullRNS-HEAAN
def ks_fullrns_heaan(sq_rho, h, sigma, k, L, Q_L, P, N):
    eta_ks_sq = (Q_L / P) * (Q_L / P) * N * sigma * sigma * (L * L + 1)
    eta_ks_sq += (k * k + 1) * (h + 1)
    eta_ks_sq *= (1.0/12)
    return sq_rho + eta_ks_sq

# Textbook CKKS rescale as implemented in HEAAN 
def rs_textbook(sq_rho, delta, h):
    return sq_rho / (delta**2) + (1.0 / 12) * (h + 1)

# Key switching for RNS variant following Lemma 10 (eprint version) in the situation of FullRNS-HEAAN
# NB: Q_ell = prod_{i=1 to q_ell} q_i for q_i the RNS moduli
def rs_fullrns_heaan(sq_rho, q_ell, h):
    return sq_rho / (q_ell**2) + (1.0 / 12) * (h + 1)

# Given a variance sq_rho, convert to an upper bound on the error coefficients with probability alpha
def rho_to_ringerror_alpha(sq_rho, N, alpha):
    return sqrt(2 * sq_rho) * erfinv(pow(1 - alpha, 1 / N))

# Given a variance sq_rho, convert to an upper bound on the total error in C^N/2
def rho_to_complexerror_alpha(sq_rho, N, alpha, delta):
    tail = sqrt(-N * log(1 - pow(1 - alpha, 2 / N)))
    return sqrt(sq_rho) * tail / delta

# Given a variance sq_rho, convert to an upper bound on the real error in C^N/2
def rho_to_realerror_alpha(sq_rho, N, alpha, delta):
    tail = erfinv(pow(1 - alpha, 2 / N))
    return sqrt(N * sq_rho) * tail / delta


# Top level function for Textbook HEAAN
def CLT_textbook(N, alpha, delta, B, sigma, q, h, tag=None):
    # For premult, need a bound on ||m||^2_2  = 2||m||^can_2 / N = Delta^2 * B^2. After add, (||m1 + m2||^can_2)^2 \leq (||m1||^can_2 + ||m2||^can_2)^2 \leq (2 * Delta * B)^2.
    message_bound = delta**2 * B**2
    var_fresh = fresh_heaan(sigma, N, h)
    var_fresh += 1. / 12  # Account for encoding noise
    var_fresh = rs_textbook(var_fresh, q, h)
    var_add = add(var_fresh, var_fresh)
    var_premult = mult(var_fresh, var_add, message_bound, 4 * message_bound, N)
    var_ks = ks_textbook(var_premult, N, h, sigma, q, q)
    var_rs = rs_textbook(var_ks, delta, h)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N, 2)
    result["log(q)"] = round(log(q, 2), 1)
    result["log(Δ)"] = log(delta, 2)
    result["α"] = alpha
    result["log(MB)"] = log(message_bound, 2)
    result["c+"] = round(log(rho_to_complexerror_alpha(var_add,
                         N, alpha, delta), 2), 2)
    result["cx"] = round(
        log(rho_to_complexerror_alpha(var_rs, N, alpha, delta), 2), 2)
    result["r+"] = round(log(rho_to_realerror_alpha(var_add,
                         N, alpha, delta), 2), 2)
    result["rx"] = round(
        log(rho_to_realerror_alpha(var_rs, N, alpha, delta), 2), 2)
    result["Rq+"] = round(log(rho_to_ringerror_alpha(var_add, N, alpha), 2), 2)
    result["Rqx"] = round(log(rho_to_ringerror_alpha(var_rs, N, alpha), 2), 2)

    return result

# Top level function \CC -> \CC for FullRNS-HEAAN
def CLT_fullrns_heaan(N, alpha, delta, B, sigma, q, k, h, L, Q_L, P, tag=None):
    message_bound = delta**2 * B**2  # B is the size bound on each entry in the message space
    var_fresh = fresh_heaan(sigma, N, h)
    var_fresh += 1. / 12  # Account for encoding noise
    var_add = add(var_fresh, var_fresh)
    var_premult = mult(var_fresh, var_add, message_bound, 4 * message_bound, N)
    var_ks = ks_fullrns_heaan(var_premult, h, sigma, k, L, Q_L, P, N)
    var_rs = rs_fullrns_heaan(var_ks, q, h)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N, 2)
    result["log(q)"] = round(log(q, 2), 1)
    result["k"] = k
    result["log(Δ)"] = log(delta, 2)
    result["α"] = alpha
    result["log(MB)"] = log(message_bound, 2)
    result["c+"] = round(log(rho_to_complexerror_alpha(var_add,
                         N, alpha, delta), 2), 2)
    result["cx"] = round(
        log(rho_to_complexerror_alpha(var_rs, N, alpha, delta), 2), 2)
    result["r+"] = round(log(rho_to_realerror_alpha(var_add,
                         N, alpha, delta), 2), 2)
    result["rx"] = round(
        log(rho_to_realerror_alpha(var_rs, N, alpha, delta), 2), 2)

    return result