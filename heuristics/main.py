"""
Structure for heuristics code files is as follows. 
- config.py: a file storing parameters
- priorCE.py: to produce prior worst-case canonical embedding bounds for CKKS as implemented in HEAAN and FullRNS-HEAAN, OpenFHE and SEAL, for comparison to our work
- CLT.py: functions needed to produce the column "CLT" in results tables for CKKS as implemented in HEAAN and FullRNS-HEAAN
- CE.py: functions needed to produce the column "CE" in results tables for CKKS as implemented in HEAAN and FullRNS-HEAAN
- WCR.py: functions needed to produce the column "WCR" in results tables for CKKS as implemented in HEAAN
- main.py: top level function
- output.txt: Output of running run_all() in main.py
"""

from config import *
from priorCE import *
from CE import *
from WCR import *
from CLT import *


def run_all(skip = {}):

    for i in range(0, len(list_of_n)):
    
        print(" -------------------------- Results for n = {} --------------------------".format(list_of_n[i]))

        if "CE" not in skip:
            print("------- CE Results -------")
            print(CE_textbook(list_of_n[i], textbook_input_bound, sigma, delta, alpha, h_heaan, list_of_Q_textbook[i], tag="CE textbook"))  #P = Q_ell
            print(CE_fullrns_heaan(list_of_n[i], sqrt(2)*input_bound, sigma, delta, alpha, h_heaan, list_of_k_fullrns_heaan[i], list_of_q_ell_fullrns_heaan[i], list_of_L_fullrns_heaan[i], list_of_Q_L_fullrns_heaan[i], list_of_P_fullrns_heaan[i], tag="CC->CC CE FullRNS-HEAAN"))
            print("\n")

        if "WCR" not in skip:
            print("------- WCR Results -------")
            print(WCR_textbook(list_of_n[i], textbook_input_bound, sigma, delta, alpha, h_heaan, list_of_Q_textbook[i], tag="WCR textbook"))  #P = Q_ell
            print("\n")

        if "CLT" not in skip:
            print("------- CLT Results -------")
            print(CLT_textbook(list_of_n[i], alpha, delta, textbook_input_bound, sigma, list_of_Q_textbook[i], h_heaan, tag="CLT textbook"))  #P = Q_ell
            print(CLT_fullrns_heaan(list_of_n[i], alpha, delta, input_bound, sigma, list_of_q_ell_fullrns_heaan[i], list_of_k_fullrns_heaan[i], h_heaan, list_of_L_fullrns_heaan[i], list_of_Q_L_fullrns_heaan[i], list_of_P_fullrns_heaan[i], tag="CLT-RNS"))
            print("\n")

        if "priorCE" not in skip:
            print("------- Prior CE Results as in [C:GHS12] etc -------")
            print(priorCE_textbook(delta, textbook_input_bound, list_of_n[i], sigma, h_heaan, list_of_Q_textbook[i], tag="priorCE textbook"))
            print(priorCE_fullrns_heaan(delta, input_bound, list_of_n[i], sigma, h_heaan, list_of_q_ell_fullrns_heaan[i], list_of_k_fullrns_heaan[i], list_of_L_fullrns_heaan[i], list_of_Q_L_fullrns_heaan[i], list_of_P_fullrns_heaan[i], tag="Prior-CE-RNS"))
        print("\n\n")
    return 0

run_all()
