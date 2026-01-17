def LB_info_set(H, p, w):
    '''
    Lee & Brickell ISD; input:
    - H: parity-check matrix
    - p: weight to be used enumeration in information set
    - w: searched weight
    '''

    #Code parameters
    n = H.ncols(); r = H.nrows(); k = n-r
        
    F2 = GF(2)

    #Sample random permutation, apply it to H
    P = Permutations(n).random_element().to_matrix().change_ring(F2)
    perm_H = H*P

    #Find invertible r x r matrix, apply permutation P_is to bring
    #it to leftmost coordinates
    I = codes.LinearCode(perm_H).information_set()
    not_I = [i for i in range(n) if i not in I]

 #   print("computing perm")
 #   print("---------- ",r,len(I),rank(H))
    #Define matrix P_is
    P_is = matrix(GF(2),n,n)
    for j in range(r):
        P_is[I[j],j] = 1
  #  print("first done")

    for j in range(k):
        P_is[not_I[j],r+j] = 1
  #  print("--------> done")

    #Apply P_is and permute again
    perm_H = perm_H*P_is
    perm_H_sys = perm_H[:,0:r]^-1 * perm_H

    #Start enumeration
  #  print("start enumeration")
    for pos in Combinations(k,p):
        #print("computing s")
        #print("pos = ",pos)

        s = sum([perm_H_sys[:,r+j] for j in pos])
        left_w = vector(s).list().count(1) #weight
        #print("------> done")

        #If weight is ok, return found codeword

        if left_w == (w-p):
            #print("desired weight found")
            c = matrix(F2,1,n) #build codeword
            for j in pos:
                c[0,r+j] = 1
            c[0,0:r] = s.transpose()

            return 1, c*(P*P_is)^-1
   # print("FAILURE")
    return 0, matrix(F2,1,n)

########################################

def sample_new_row(H, u, v, p):
    '''
    Samples vector c weight v such that H[0:u, :]*c^T = 0
    u is the number of rows we need to consider from H
    p is the parameter for Lee and Brickell ISD
    '''
    n = H.ncols()
    C = codes.LinearCode(H[0:u, :])
    H_isd = matrix(GF(2), u+1, n)
    H_isd[0:u, :] = H[0:u,:]

    for i in range(n):
        H_isd[u, i] = 1

    #Recompute H to handle cases where the all ones vector is already in the code
    H_isd = codes.LinearCode(H_isd).generator_matrix()

    #start looping until a valid new row is found
    row_found = 0
    num_iter = 0
    while row_found == 0:

        #call ISD until a codeword with the desired weight is found
        ok_isd = 0
        while ok_isd == 0:
            if ((v-0.5*(u+1))>p):
                ok_isd, c = LB_info_set(H_isd, v-ceil(0.5*(u+1)), v)
                num_iter += 1
            else:
                ok_isd, c = LB_info_set(H_isd, p, v)
                num_iter += 1

        if ok_isd:
            H[u,:] = c
           # print("-------->",rank(H))
            if rank(H[0:u+1,:]) == (u+1):
                row_found = 1

             #   print("----------------> RANK DEFICIENCY")
        if num_iter == 100:
            return num_iter, H

    return num_iter, H
