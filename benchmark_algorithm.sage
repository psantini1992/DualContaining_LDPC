reset()
load('LB.sage')
load('utils.sage')

n = 150
k = 80
v = 10

d_rnd = gv_rnd(n, k)
print("n = ",n,", k = ",k,", GV distance = ",d_rnd)
ewd = exact_distrib(n,n-k-1,v)

print("LDPC: v = ",v,", m_v = ",ewd[v])

p = 3 #parameter for LB ISD

num_test = 100

import time
import datetime

r = n-k #code redundancy
F2 = GF(2)

avg_time = 0
avg_num_iter = 0

success_rate = 0
dimensions = vector(ZZ,r+1)

for idx in range(num_test):

    success = 1

    start = time.time()

    #create matrix
    H = matrix(F2,r,n)
    pos = Combinations(n,v).random_element()
    for j in pos:
        H[0,j] = 1

    #Sample another row of weight v,
    #until it is orthogonal to the already found row
    ok = 0
    while ok == 0:
        pos = Combinations(n,v).random_element()
        x = sum([H[0,j] for j in pos])
        if x == 0:
            ok = 1
            for j in pos:
                H[1,j] = 1

    this_r = 2
    avg_num_iter += 2

    print("Code #",idx,", simulation started at ",datetime.datetime.now())

    while this_r < r:
        num_iter, H = sample_new_row(H, this_r, v, p)
        if num_iter == 100:
            dimensions[this_r] += 1
            print("\r\nFAILURE, TOO MANY ITERATIONS!")
            this_r = r
            success = 0
        avg_num_iter += num_iter
        this_r += 1
        sys.stdout.write("\r--> Done: u = %i" % this_r)
        sys.stdout.flush()

    success_rate += success

    end = time.time()
    print("\r\n--> Done: time = ",end-start,"s")
    print("--> Avg num iter = ",N(avg_num_iter/(r*(1+idx))))
    print("============================================")
    print(" ")
    avg_time += (end-start)

print("SIMULATION DONE")
print("- n = ",n)
print("- r = ",r)
print("- v = ",v)
print("- m_v = ",ewd[v])
print("- GV distance = ",d_rnd)
print("- LB with p = ",p)
print("- num test = ",num_test)
print("- avg time = ",N(avg_time/num_test))
print("- avg num iter = ",N(avg_num_iter/(r*num_test)))
print("- success rate = ",N(success_rate/(num_test)))
