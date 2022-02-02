######################
# CLT CKKS ESTIMATES #
######################
from math import sqrt, log, floor
from scipy.special import erfinv,erf
from config import *
# the CLT approach tracks the growth of the variance of the error coefficients in the ring
# does not generate bounds directly -- once a final rho is calculated, translate to worst
# case bound on error coeffs, total error in C^N/2, or real error in C^N/2
# require independence between errors to be valid.
# each function returns a value for rho^2 i.e. the variance


def clt_fresh(sigma, N):
    return sigma * sigma * (1 + 4 * N / 3)


def clt_add(sq_rho_1, sq_rho_2):  # need independence
    return sq_rho_1 + sq_rho_2


def clt_mult(sq_rho_1, sq_rho_2, M_1, M_2, N):
    # M_1 and M_2 are (bounds on) |m_1|^2, |m_2|^2 the squared 2 norm of the polynomial
    return N * sq_rho_1 * sq_rho_2 + sq_rho_1 * M_2 + sq_rho_2 * M_1


def clt_ks(sq_rho, N, sigma):
    # removing the q_l from the formula bc gets deleted below anyway
    return sq_rho + sigma**2 * 3 * N


def clt_rs(sq_rho, delta, N):
    # rescale a random v with variance sq_rho by delta
    return sq_rho / (delta**2) + 1 / 12 + N / 18


def rho_to_ringerror(sq_rho):
    # given a variance sq_rho, convert to an upper bound on the error coefficients with probability alpha
    return 6 * sqrt(sq_rho)


def rho_to_ringerror_alpha(sq_rho, N, alpha):
    # given a variance sq_rho, convert to an upper bound on the error coefficients with probability alpha
    return sqrt(2 * sq_rho) * erfinv(pow(1 - alpha, 1 / N))


def rho_to_complexerror_alpha(sq_rho, N, alpha, delta):
    # translate to total error in complex
    tail = sqrt(-N * log(1 - pow(1 - alpha, 2 / N))) 
    return sqrt(sq_rho) * tail / delta


def rho_to_realerror_alpha(sq_rho, N, alpha, delta):
    # a bound on just the real error
    tail = erfinv(pow(1 - alpha, 2 / N))
    return sqrt(N * sq_rho) * tail / delta


# We go through the following circuit: generate two
# fresh ciphertexts, add them, multiply a fresh ct
# with the result of the addition, modulus switch.

def test_clt_bounds(N, alpha, delta, B, sigma, q, tag=None):
    message_bound = N * B * B
    var_fresh = clt_fresh(sigma, N)
    var_add = clt_add(var_fresh, var_fresh)
    var_mult = clt_mult(var_fresh, var_add, message_bound, 4 * message_bound, N)
    var_mult = clt_ks(var_mult, N, sigma)
    var_rs = clt_rs(var_mult, delta, N)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N, 2)
    result["log(q)"] = round(log(q, 2), 1)
    result["log(Δ)"] = log(delta, 2)
    result["α"] = None
    result["fresh noise"] = round(log(rho_to_ringerror(var_fresh), 2), 2)
    result["log(MB)"] = round(log(message_bound, 2), 1)
    result["+"] = round(log(rho_to_ringerror(var_add), 2), 2)
    result["x"] = round(log(rho_to_ringerror(var_rs), 2), 2)

    return result


def test_clt_bounds_alpha(N, alpha, delta, B, sigma, q, tag=None):
    message_bound = N * B * B
    var_fresh = clt_fresh(sigma, N)
    var_add = clt_add(var_fresh, var_fresh)
    var_mult = clt_mult(var_fresh, var_add, message_bound, 4 * message_bound, N)
    var_mult = clt_ks(var_mult, N, sigma)
    var_rs = clt_rs(var_mult, delta, N)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N, 2)
    result["log(q)"] = round(log(q, 2), 1)
    result["log(Δ)"] = log(delta, 2)
    result["α"] = alpha
    result["log(MB)"] = round(log(message_bound, 2), 1)
    result["+"] = round(log(rho_to_ringerror_alpha(var_add, N, alpha), 2), 2)
    result["x"] = round(log(rho_to_ringerror_alpha(var_rs, N, alpha), 2), 2)

    return result

# now we do the same for encoding

def clt_encode_message_bound(input_bound, delta):
    return delta**2 * input_bound**2

# this is a duplicate of rho to complex error alpha, line 49
def clt_decode_complex(rho, N, alpha, delta):
    blob = sqrt(- N * log(1 - alpha**(2 / N), 2))
    return rho * blob / delta

## this is a duplicate of rho to real error alpha, line 55
def clt_decode_real(rho, N, alpha, delta):
    return sqrt(rho * N) * erfinv((1 - alpha)**(2 / N)) / delta


def test_clt_bounds_decoding_alpha(N, alpha, delta, B, sigma, q, tag=None):
    message_bound = clt_encode_message_bound(B, delta)
    var_fresh = clt_fresh(sigma, N)
    var_fresh += 1. / 12
    var_add = clt_add(var_fresh, var_fresh)
    var_mult = clt_mult(var_fresh, var_add, message_bound, 4 * message_bound, N)
    var_mult = clt_ks(var_mult, N, sigma)
    var_rs = clt_rs(var_mult, delta, N)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N, 2)
    result["log(q)"] = round(log(q, 2), 1)
    result["log(Δ)"] = log(delta, 2)
    result["α"] = alpha
    result["log(MB)"] = log(message_bound, 2)
    result["c+"] = round(log(rho_to_complexerror_alpha(var_add, N, alpha, delta), 2), 2)
    result["cx"] = round(log(rho_to_complexerror_alpha(var_rs, N, alpha, delta), 2), 2)
    result["r+"] = round(log(rho_to_realerror_alpha(var_add, N, alpha, delta), 2), 2)
    result["rx"] = round(log(rho_to_realerror_alpha(var_rs, N, alpha, delta), 2), 2)

    return result
