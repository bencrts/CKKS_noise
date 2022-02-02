######################################
# CANONICAL EMBEDDING CKKS ESTIMATES #
######################################
from math import *
from config import *

# This script will heuristically estimate the following homomorphic
# circuits:
#ct4 = ct1 + ct2
#ct5 = ct3 * ct4
#ct6 = Rescale(ct5)
# where ct1, ct2, and ct3 are fresh ciphertexts. We then compare this against
# real noise growth using the HEAAN library

# Noise growth for individual operations

def ckks_fresh(sigma, N):
    return 6 * sigma * sqrt(N * (4 * N / 3 + 1))


def ckks_fresh_alpha(sigma, N, alpha):
    nu = log(1 - pow(1 - alpha, 2. / N))
    return sigma * sqrt(N * (4 * N / 3 + 1)) * sqrt(-nu)


def ckks_add(input_noise1, input_noise2):
    return input_noise1 + input_noise2


def ckks_ks(N, sigma):
    return sigma * sqrt(3) * N


def ckks_ks_alpha(N, sigma, alpha):
    nu = log(1 - pow(1 - alpha, 2 / N))
    return sigma * N * sqrt(-nu) * sqrt(1 / 12)


def ckks_scale(N):
    return 6 * sqrt(N / 12 + N**2 / 18)


def ckks_scale_alpha(N, alpha):
    nu = log(1 - pow(1 - alpha, 2 / N))
    return sqrt(N / 12 + N**2 / 18) * sqrt(-nu)

# ckks_mult includes key switch but not scale

def ckks_mult(input_noise1, input_noise2, message_bound1, message_bound2, N, sigma, delta):
    # blob = (1/P)*q_current*ckks_ks(n, q_current, sigma) old
    blob = ckks_ks(N, sigma)
    return (input_noise1 * message_bound2 + input_noise2 * message_bound1 +
            input_noise1 * input_noise2 + blob) / delta + ckks_scale(N)


def ckks_mult_alpha(input_noise1, input_noise2, message_bound1, message_bound2, N, sigma, alpha, delta):
    blob = ckks_ks_alpha(N, sigma, alpha)
    return (input_noise1 * message_bound2 + input_noise2 * message_bound1 +
            input_noise1 * input_noise2 + blob) / delta + ckks_scale_alpha(N, alpha)

# Now let's get the actual noise growth in in a circuit

def get_can_emb_noise(input_noise, q_current):
    noise_mod_q = math.fmod(input_noise, q_current)
    return noise_mod_q

# We go through the following circuit: generate two
# fresh ciphertexts, add them, multiply a fresh ct
# with the result of the addition, modulus switch.

def test_can_emb_noise(N, B, sigma, delta, tag = None):
    # This is the 6rootV version of the canonical embedding, ring to ring, circuit
    message_bound = 2*B*N/pi
    fresh = ckks_fresh(sigma, N)
    noise_add = ckks_add(fresh, fresh)
    noise_mult = ckks_mult(noise_add, fresh, 2 * message_bound, message_bound, N, sigma, delta)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N,2)
    result["log(Δ)"] = log(delta, 2)
    result["α"] = alpha
    result["log(MB)"] = log(message_bound, 2)
    result["fresh-noise"] = round(log(fresh, 2), 1)
    result["+"] = round(log(noise_add, 2), 1)
    result["x"] = round(log(noise_mult, 2), 1)

    return result


def test_can_emb_noise_alpha(N, B, sigma, delta, alpha, tag = None):
    message_bound = 2 * B * N / pi
    fresh = ckks_fresh_alpha(sigma, N, alpha)
    noise_add = ckks_add(fresh, fresh)
    noise_mult = ckks_mult_alpha(noise_add, fresh, 2 * message_bound, message_bound, N, sigma, alpha, delta)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N,2)
    result["log(Δ)"] = log(delta, 2)
    result["α"] = alpha
    result["log(MB)"] = round(log(message_bound, 2),2)
    result["+"] = round(log(noise_add, 2), 2)
    result["x"] = round(log(noise_mult, 2), 2)

    return result

# now we do the same for encoding

def can_emb_encode(N):
    return N / pi


def can_emb_decode(noise, delta):
    return noise / delta


def test_can_emb_encoding_alpha(N, B, sigma, delta, alpha, tag = None):
    message_bound = B * delta
    fresh = ckks_fresh_alpha(sigma, N, alpha)
    fresh += can_emb_encode(N)
    noise_add = ckks_add(fresh, fresh)
    noise_mult = ckks_mult_alpha(noise_add, fresh, 2* message_bound, message_bound, N, sigma, alpha, delta)
    decode_add = can_emb_decode(noise_add, delta)
    decode_mult = can_emb_decode(noise_mult, delta)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N,2)
    result["log(Δ)"] = log(delta, 2)
    result["α"] = alpha
    result["log(MB)"] = log(message_bound, 2)
    result["+"] = round(log(decode_add, 2), 2)
    result["x"] = round(log(decode_mult, 2), 2)

    return result

# going to write a version of above, which uses 6rootV estimates instead of alpha

def test_can_emb_encoding(N, input_bound, sigma, delta, tag = None):
    message_bound = input_bound * delta
    fresh = ckks_fresh(sigma, N)
    fresh += can_emb_encode(N)
    fresh_prec = can_emb_decode(fresh,delta)    
    noise_add = ckks_add(fresh, fresh)
    noise_mult = ckks_mult(noise_add, fresh, 2 * message_bound, message_bound, N, sigma, delta)
    decode_add = can_emb_decode(noise_add, delta)
    decode_mult = can_emb_decode(noise_mult, delta)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N,2)
    result["log(Δ)"] = log(delta, 2)
    result["fresh-prec-loss"] = round(log(fresh_prec,2),2)
    result["α"] = alpha
    result["log(MB)"] = log(message_bound, 2)
    result["+"] = round(log(decode_add, 2), 2)
    result["x"] = round(log(decode_mult, 2), 2)

    return result
