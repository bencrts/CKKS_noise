from math import sqrt

# Value of delta used in all experiments
delta = 2**40

# Value of alpha used in all experiments
alpha = 0.00001

# Value of input bound (in the complex space) for all experiments
# The HEAAN and FullRNS-HEAAN inputs are chosen to have both real and imaginary components in [0, 1]
textbook_input_bound = sqrt(2)
input_bound = sqrt(2)

# Values of n for all experiments
n_4096 = 4096  # 12 bits
n_8192 = 8192  # 13 bits
n_16384 = 16384  # 14 bits
n_32768 = 32768  # 15 bits
list_of_n = [n_4096, n_8192, n_16384, n_32768]

# Value of sigma for all experiments
sigma = 3.2

# Values of Q_ell used in textbook CKKS experiments (HEAAN)
Q_4096_textbook = 18014398509481984  # 54 bits
Q_8192_textbook = 649037107316853453566312041152512  # 109 bits
Q_16384_textbook = 842498333348457493583344221469363458551160763204392890034487820288  # 219 bits
Q_32768_textbook = 22713710134237715329666368996500141698551292521478689383796568724394977753543685103943470334805111423773828800195818060422956300894208  # 443 bits
list_of_Q_textbook = [Q_4096_textbook, Q_8192_textbook, Q_16384_textbook, Q_32768_textbook]

# Values of q_ell used in FullRNS-HEAAN
list_of_q_ell_fullrns_heaan = [1099511922689, 1099511480321, 1099514314753, 1099516870657]

# Values of Q_ell used in FullRNS-HEAAN
Q_L_4096_RNS = 2535300860448494974062284464129
Q_L_8192_RNS = 2535300860448494974062284464129
Q_L_16384_RNS = 3370001669392531206832073295534520426823353513034166032940060147713
Q_L_32768_RNS = 5415430856483162131116411537630132849886909475854026828679402461173808845823023195972815525471370082192395894955789727182487553
list_of_Q_L_fullrns_heaan = [Q_L_4096_RNS, Q_L_8192_RNS, Q_L_16384_RNS, Q_L_32768_RNS]

# Value of h used in textbook CKKS experiments (HEAAN) and used in FullRNS-HEAAN
h_heaan = 64

# Values of k used in FullRNS-HEAAN
list_of_k_fullrns_heaan = [3, 3, 6, 11]

# Values of L used in FullRNS-HEAAN
list_of_L_fullrns_heaan = [2, 2, 5, 10]

# Values of P used in FUllRNS-HEAAN
list_of_P_fullrns_heaan = [1329228273087045263862449605477367809,
1329229580352387071335847315304087553,
1766859912869502984729736217977187591167433174269945269833750071089954817,
2839147597473542536594182342178826719455684754759561317787248478193202993464905588804670881347564059435619411970833615385437813669889]
