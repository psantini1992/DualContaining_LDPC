load('LB.sage')
load('utils.sage')

n = 80 #code length
r = 40 #code redundancy
v = 7 #row weight

F2 = GF(2)
k = n-r #code dimension

d = gv_rnd(n, k)
print("GV distance: d = ",d)

m_w_coeffs = vector(RR, d+1)
for w in range(1,d+1):
    pr = sum([1.*binomial(v,x)*binomial(n-v,w-x)/binomial(n,w) for x in range(0,w+1,2)])
    m_w = binomial(n,w)*pr^(n-k)
    print("w = ",w,", m_w = ", m_w)
    m_w_coeffs[w] = m_w

w_min = 4 #min weight
w_max = 10 #max weight
num_test = 200 #number of ISD calls
num_codes = 200 #number of considered codes
p = 3 #parameter for Lee&Brickell

#Values of m_w
th_m_w_filename = "Results/test_ISD/th_m_w_n"+str(n)+"_r"+str(r)+"_v"+str(v)+"_p"+str(p)+".txt"
emp_m_w_filename = "Results/test_ISD/emp_m_w_n"+str(n)+"_r"+str(r)+"_v"+str(v)+"_p"+str(p)+".txt"

#Number of iterations
th_num_iteration_filename = "Results/test_ISD/th_num_iteration_n"+str(n)+"_r"+str(r)+"_v"+str(v)+"_p"+str(p)+".txt"
emp_num_iteration_filename = "Results/test_ISD/emp_num_iteration_n"+str(n)+"_r"+str(r)+"_v"+str(v)+"_p"+str(p)+".txt"

#Number of codes with codewords of weight w
emp_num_codes_filename = "Results/test_ISD/emp_num_codes_n"+str(n)+"_r"+str(r)+"_v"+str(v)+"_p"+str(p)+".txt"


with open(th_m_w_filename,'w') as f:
    f.write("X Y\r\n")
with open(emp_m_w_filename,'w') as f:
    f.write("X Y\r\n")
with open(th_num_iteration_filename,'w') as f:
    f.write("X Y\r\n")
with open(emp_num_iteration_filename,'w') as f:
    f.write("X Y\r\n")
with open(emp_num_codes_filename,'w') as f:
    f.write("X Y\r\n")

#Values of m_w
th_m_w_vals = []
emp_m_w_vals = []

#Number of iterations
th_num_iteration = []
emp_num_iteration = []

#Number of codes with codewords of weight w
num_codes_with_w = []

for w in range(w_min, w_max+1):


    th_success_pr = min(1, N(m_w_coeffs[w]*binomial(k,p)*binomial(n-k,w-p)/binomial(n,w)))
    print("w = ",w,", Exp. Num Iter = ", 1/th_success_pr)

    th_num_iteration.append((w, 1/th_success_pr))
    th_m_w_vals.append((w, m_w_coeffs[w]))

    #Empirical values
    num_iter = 0; num_codes_ok = 0; emp_m_w = 0; num_success = 0

    for num in range(num_codes):

        H = sample_H(n, r, v) #sample parity-check matrix
        found_c = [] #found codewords

        #Call ISD for num_test times
        for num_isd in range(num_test):
            ok, c = LB_info_set(H, p, w)

            #if new codeword is found, append it to found codewords
            if ok:
                num_success += 1

                if (c not in found_c):
                    found_c.append(c)
            if (num_isd%20) == 0:
                print("--> w = ",w,", Code #",num+1,", "+str(round((num_isd+1)/num_test*100,2))+"% Done: Th m_w = ",N(m_w_coeffs[w]),", Num found codewords = ",len(found_c))

        if len(found_c)>0:
            num_codes_ok += 1

        emp_m_w += len(found_c)
        print("w = ",w,", Code #",num+1," Completed: Th m_w = ",N(m_w_coeffs[w]),", Emp. m_w = ",N(emp_m_w/(num+1)))

    #update values
    num_codes_with_w.append((w, num_codes_ok))
    emp_m_w_vals.append((w, N(emp_m_w/(num_codes))))
    emp_num_iteration.append((w, N(num_codes*num_test/num_success)))

    with open(th_m_w_filename,'a') as f:
        f.write(str(w)+" "+str(m_w_coeffs[w])+"\r\n")
    with open(emp_m_w_filename,'a') as f:
        f.write(str(w)+" "+str(N(emp_m_w/num_codes))+"\r\n")
    with open(th_num_iteration_filename,'a') as f:
        f.write(str(w)+" "+str(1/th_success_pr)+"\r\n")
    with open(emp_num_iteration_filename,'a') as f:
        f.write(str(w)+" "+str(N(num_codes*num_test/num_success))+"\r\n")
    with open(emp_num_codes_filename,'a') as f:
        f.write(str(w)+" "+str(N(num_codes_ok/num_codes))+"\r\n")
