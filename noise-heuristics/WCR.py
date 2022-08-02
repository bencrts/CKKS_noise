#########################################
# WORST CASE IN THE RING CKKS ESTIMATES #
#########################################
from math import *
from scipy.special import erfinv
from config import *
# The set-up is as follows: this script will
# heuristically estimate the following homomorphic
# circuit: ct1, ct2, ct3 are fresh. We do
#ct4 = ct1 + ct2
#ct5 = ct3*ct4
#ct6 = Rescale(ct5)
# Then we compare this against real noise growth in
#experiments in HEEAN

# Note to self: in the previous experiements we had
# to multiply all noises by q_current for this to work.
# Probably will need to do this again but let's see.

# Let's get the functions that show how much noise growth
# each homomorphic operation incurs


def ckks_fresh(sigma, N):
    return 6 * sigma * sqrt(1 + (4 * N) / 3)


def ckks_fresh_alpha(sigma, N, alpha):
    return sigma * sqrt(2 + (8 * N) / 3) * erfinv(pow(1 - alpha, 1 / N))


def ckks_add(input_noise1, input_noise2):
    return input_noise1 + input_noise2


def ckks_ks(N, sigma):
    return sigma * sqrt(3 * N)


def ckks_ks_alpha(N, sigma, alpha):
    return sigma * sqrt(N / 6) * erfinv(pow(1 - alpha, 1 / N))


def ckks_scale(N):
    return sqrt(3 + 2 * N)


def ckks_scale_alpha(N, alpha):
    return sqrt((3 + 2 * N) / 18) * erfinv(pow(1 - alpha, 1 / N))

# ckks_mult includes key switch but not scale


def ckks_mult(input_noise1, input_noise2, message_bound1, message_bound2, N, sigma, delta):
    # blob = (1/P)*q_current*ckks_ks(n, q_current, sigma) old
    blob = ckks_ks(N, sigma)

    # P^{-1} * q_l = 1 in our exps
    return (N * input_noise1 * message_bound2 + N * input_noise2 * message_bound1 +
            N * input_noise1 * input_noise2 + blob) / delta + ckks_scale(N)


def ckks_mult_alpha(input_noise1, input_noise2, message_bound1, message_bound2, N, sigma, alpha, delta):
    blob = ckks_ks_alpha(N, sigma, alpha)
    return N * (input_noise1 * message_bound2 + input_noise2 * message_bound1 +
                input_noise1 * input_noise2) / delta + blob + ckks_scale_alpha(N, alpha)

# Now let's get the actual noise growth in in a circuit (actual size, not bits)


def get_lwe_noise(input_noise, q_current):
    noise_mod_q = input_noise % q_current
    return noise_mod_q


# We go through the following circuit: generate two
# fresh ciphertexts, add them, multiply a fresh ct
# with the result of the addition, modulus switch.

def test_lwe_noise(N, B, sigma, delta, tag=None):
    message_bound = B
    fresh = ckks_fresh(sigma, N)
    noise_add_C1 = ckks_add(fresh, fresh)
    noise_add_mult_C3 = ckks_mult(noise_add_C1, fresh, 2 * message_bound, message_bound, N, sigma, delta)
    result = {}
    result["type"] = tag
    result["log(N)"] = log(N, 2)
    result["log(Δ)"] = log(delta, 2)
    result["α"] = None
    result["fresh noise"] = round(log(fresh, 2), 2)
    result["log(MB)"] = round(log(message_bound, 2), 2)
    result["+"] = round(log(noise_add_C1, 2), 2)
    result["x"] = round(log(noise_add_mult_C3, 2), 2)
    return result


def test_lwe_noise_alpha(N, B, sigma, delta, alpha, tag=None):
    message_bound = B
    fresh = ckks_fresh_alpha(sigma, N, alpha)
    noise_add_C1 = ckks_add(fresh, fresh)
    noise_add_mult_C3 = ckks_mult_alpha(noise_add_C1, fresh, 2 * message_bound, message_bound, N, sigma, alpha, delta)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N, 2)
    result["log(Δ)"] = log(delta, 2)
    result["α"] = alpha
    result["log(MB)"] = round(log(message_bound, 2), 2)
    result["+"] = round(log(noise_add_C1, 2), 2)
    result["x"] = round(log(noise_add_mult_C3, 2), 2)

    return result

# now we do the same for encoding


def wcr_message_bound(B, delta):
    return B * delta


def wcr_decode(N, noise, delta):
    return (2 * N * noise) / (pi * delta)


def test_lwe_noise_encoding_alpha(N, B, sigma, delta, alpha, tag=None):
    message_bound = wcr_message_bound(B, delta)
    fresh = ckks_fresh_alpha(sigma, N, alpha)
    fresh += 0.5
    noise_add_C1 = ckks_add(fresh, fresh)
    noise_add_mult_C3 = ckks_mult_alpha(noise_add_C1, fresh, 2 * message_bound, message_bound, N, sigma, alpha, delta)
    decode_add_C1 = wcr_decode(N, noise_add_C1, delta)
    decode_add_mult_C3 = wcr_decode(N, noise_add_mult_C3, delta)

    result = {}
    result["type"] = tag
    result["log(N)"] = log(N, 2)
    result["log(Δ)"] = log(delta, 2)
    result["α"] = alpha
    result["log(MB)"] = round(log(message_bound, 2), 2)
    result["+"] = round(log(decode_add_C1, 2), 2)
    result["x"] = round(log(decode_add_mult_C3, 2), 2)

    return result
