def sample_H(n, r, v):
    '''
    Sample random matrix with n columns, r rows and each row of weight v
    '''
    
    #sample H
    H = matrix(GF(2), r, n)
    for i in range(r):
        pos = Combinations(n,v).random_element()
        for j in pos:
            H[i,j] = 1
    
    return H

def exact_distrib(n,r,v):
    '''
    Weight distribution for codes with parity-check matrix wih constant row weight v
    '''
    coeffs = vector(RR,n+1)
    for w in range(1,n+1):
        pr = sum([N(binomial(v,i)*binomial(n-v,w-i)/binomial(n,w)) for i in range(0,min(v,w)+1,2)])
        m_w = log(1.*binomial(n,w),2)+r*log(pr,2)
        coeffs[w] = 2**m_w
    
    return coeffs

def gv_rnd(n, k):
    w = 1
    while (binomial(n,w)*2^(k-n))<1:
        w += 1
    return w

def LB_cost(n,k,w,num_w):
    best_cost = 1000000000000
    best_p = 0
    for p in range(1,min(10,w)):
        cost = log(1.*(n-k)^2*(n+k)+n*binomial(k,p),2)+log(1.*binomial(n,w)/(binomial(k,p)*binomial(n-k,w-p)),2)-log(num_w,2)
        if cost < best_cost:
            best_cost = cost
            best_p = p
    return best_p, best_cost
