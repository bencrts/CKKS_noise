from WCR import *
from CLT import *
from CE import * 
from config import *


def run_all(skip = {}):

    for i in range(1, len(list_of_n)):

        if "CE" not in skip:
            print("------- CE Results -------")
            print(test_can_emb_noise(list_of_n[i], delta, sigma, delta, tag = "Rq -> Rq CE 6√V"))
            print(test_can_emb_noise_alpha(list_of_n[i], delta, sigma, delta, alpha, tag = "Rq -> Rq CE"))
            # for real to real, the complex entries have input bound 1
            print(test_can_emb_encoding_alpha(list_of_n[i], input_bound, sigma, delta, alpha, tag = "RR->RR CE"))
            # for complex to complex, the complex entries have input bound sqrt(2)
            print(test_can_emb_encoding_alpha(list_of_n[i], sqrt(2)*input_bound, sigma, delta, alpha, tag = "CC->CC CE"))
            print(test_can_emb_encoding(list_of_n[i], sqrt(2)*input_bound, sigma, delta, tag = "CC->CC CE 6√V"))
            print("\n\n")

        if "WCR" not in skip:
            print("------- WCR Results -------")
            print(test_lwe_noise(list_of_n[i],delta,sigma,delta,tag="Rq->Rq WCR 6rootV"))
            print(test_lwe_noise_alpha(list_of_n[i], delta, sigma, delta, alpha, tag = "Rq->Rq WCR"))
            print(test_lwe_noise_encoding_alpha(list_of_n[i], input_bound, sigma, delta, alpha, tag = "RR->RR WCR"))
            # for complex to complex, the complex entries have input bound sqrt(2)
            print(test_lwe_noise_encoding_alpha(list_of_n[i], sqrt(2)*input_bound, sigma, delta, alpha, tag = "CC->CC WCR"))
            print("\n\n")

        if "CLT" not in skip:
            print("------- CLT Results -------")
            print(test_clt_bounds(list_of_n[i], alpha, delta, B, sigma, list_of_q[i], tag="No-α CLT"))
            print(test_clt_bounds_alpha(list_of_n[i], alpha, delta, B, sigma, list_of_q[i], tag="Rq->Rq CLT"))
            # for real to real, the complex entries have input bound 1
            print(test_clt_bounds_decoding_alpha(list_of_n[i], alpha, delta, input_bound, sigma, list_of_q[i], tag="RR->RR CLT"))
            # for complex to complex, the complex entries have input bound sqrt(2)
            print(test_clt_bounds_decoding_alpha(list_of_n[i], alpha, delta, input_bound * sqrt(2), sigma, list_of_q[i], tag = "CC->CC CLT"))

        return 0

run_all()