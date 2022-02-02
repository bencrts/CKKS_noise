def m_k(k, a, b, Delta, N):

    T_1 = -8/(a**2 + 6*a*b + b**2)
    T_0 = 8*(a+b)/(a**2 + 6*a*b + b**2)

    if k == 0:
        return Delta * (T_0 - T_1*a)#a-bX is decreasing in X
    
    else:
        return Delta/a #2X - dX^2 has maximum at 1/d < 1/a


def Pinv_q_l(l, logDelta):
    # treating l as a count of the number of rescales done i.e. levels are increasing with depth
    # only want P^-1 q_l and for HEAAN P=2^logQ so
    return 2**(-l*logDelta)


def logA(r, a, b):

    t = (b-a)**2 / (a**2 + 6*a*b + b**2)

    return 2**r  * log(t) - log(a)


def Br_WCR(r, a, b, Delta, N, B_fresh, B_ks, B_rs):
    #ring bound 
    if r == 0:
        T_1 = -8/(a**2 + 6*a*b + b**2)
        T_0 = 8*(a+b)/(a**2 + 6*a*b + b**2)
        return b * abs(round(Delta * T_1) - Delta * T_1) + (1/Delta) * round(Delta * T_1) * B_fresh + abs(round(Delta * T_0) - Delta * T_0)+B_rs

    else:
        mk_oo = m_k(r-1, a, b, Delta, N)
        B_prev = Br_WCR(r - 1, a, b, Delta, N, B_fresh, B_ks, B_rs)
        # if on iteration r, must have performed 1 + 2(r-1) rescales so
        ks_coeff = Pinv_q_l(2*r -1, Delta)
        B_square = (1/Delta) * (N * B_prev * (B_prev + 2 * mk_oo) + ks_coeff * B_ks) + B_rs
        term1 = 2 * B_prev
        #we're now one level up, so change ks_coeff
        ks_coeff = Pinv_q_l(2*r,Delta)
        term2 = (1/Delta) * (N * (B_fresh * B_square + Delta*b* B_square + (1/Delta)*(mk_oo**2)*B_fresh) + ks_coeff * B_ks) + B_rs
        return term1 + term2


def Beta_r_WCR(r, a, b, Delta, N, B_fresh, B_ks, B_rs):
    #log noise in C^N/2 for WCR
    Bk = Br_WCR(r, a, b, Delta, N, e_fresh, e_ks, e_rs)

    return log(2 * N * Bk) - log(pi * Delta)

def Br_CE(r, a, b, Delta, N, B_fresh, B_ks, B_rs):
    #canonical embedding bound
    if r == 0:
        T_1 = -8/(a**2 + 6*a*b + b**2)
        T_0 = 8*(a+b)/(a**2 + 6*a*b + b**2)
        return b * abs(round(Delta * T_1) - Delta * T_1) + (1/Delta) * round(Delta * T_1) * B_fresh + abs(round(Delta * T_0) - Delta * T_0)+B_rs

    else:
        #need norm of previous poly
        mk_oo = m_k(r-1, a, b, Delta, N)
        B_prev = Br_CE(r - 1, a, b, Delta, N, B_fresh, B_ks, B_rs)
        # if on iteration r, must have performed 1 + 2(r-1) rescales so
        ks_coeff = Pinv_q_l(2*r -1, Delta)
        B_square = (1/Delta) * (B_prev * (B_prev + 2 * mk_oo) + ks_coeff * B_ks) + B_rs
        term1 = 2 * B_prev
        #we're now one level up, so change ks_coeff
        ks_coeff = Pinv_q_l(2*r,Delta)
        term2 = (1/Delta) * ((B_fresh * B_square + Delta*b* B_square + (1/Delta)*(mk_oo**2)*B_fresh) + ks_coeff * B_ks) + B_rs
        return term1 + term2

def Beta_r_CE(r, a, b, Delta, N, B_fresh, B_ks, B_rs):
    Bk = Br_CE(r, a, b, Delta, N, e_fresh, e_ks, e_rs)
    return log(Bk) - log(Delta)
def testit():

    sigma = 3.2
    q = 2**120

    kwds = {        
        "a": 1,
        "b": 5,
        "Delta": 2**25,
        "N": 2**12,
        "B_fresh": 6 * sigma * sqrt(4 * 2**12 /3 + 1),
        # TODO add ks properly: think this is now done?
        "B_ks": sqrt(3 + 2 * 2**12),
        "B_rs": sqrt(3 + 2 * 2**12),
    } 

    kwds2 = {
        "a": 1,
        "b": 5,
    }

    loga = []
    logbeta = []

    for r in range(10):
        loga.append(logA(r = r, **kwds2))
        logbeta.append(Beta_r_WCR(r = r, **kwds).n())
    
    return loga, logbeta

#for now doing 6rootV
def WCR_fresh(N,sigma):
    return 6*sigma*sqrt(1 + (4*N)/3)

def WCR_ks(N, sigma):
    return sigma*sqrt(3*N)

def WCR_rs(N):
    return sqrt(3+2*N)

def CE_fresh(N,sigma):
    return 6*sqrt(sigma**2 * N * (4*N/3 + 1))

def CE_ks(N,sigma):
    return sigma*sqrt(3)*N

def CE_rs(N):
    return sqrt(3*N + 2*N**2)

def identify_CP(a,b,N,Delta,sigma,loop):

    loga = []
    for r in range(loop):
        loga.append(logA(r = r, a,b))

    #get the WCR bounds
    B_fresh = WCR_fresh(N,sigma)
    B_ks = WCR_ks(N,sigma)
    B_rs = WCR_rs(N)

    for r in range(loop):
        beta = Beta_r_WCR(r, a, b, Delta, N, B_fresh, B_ks, B_rs)
        if(beta > loga[r]):
            print("WCR critical point reached at iteration" r "\n")
            exit for
    
    #get the CE bounds
    B_fresh = CE_fresh(N,sigma)
    B_ks = CE_ks(N,sigma)
    B_rs = CE_rs(N)

    for r in range(loop):
        beta = Beta_r_CE(r, a, b, Delta, N, B_fresh, B_ks, B_rs)
        if(beta > loga[r]):
            print("CE critical point reached at iteration" r "\n")
            exit for






