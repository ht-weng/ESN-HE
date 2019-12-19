import numpy as np

a = np.matrix([[1, 2, 3], 
              [4, 5, 6], 
              [7, 8, 9]])
b = np.matrix([[11, 12, 13], 
              [14, 15, 16], 
              [17, 18, 19]])

a_b = a.dot(b)

print("True Result:")
print(a_b)

a_text = a.A1
b_text = b.A1

def rot(ct, l):
    return np.roll(ct, l)

def add(a, b):
    return a + b

def mult(a, b):
    return a * b

def cmult(a, b):
    return a * b

def genUs(d):
    mat = np.zeros((d*d, d*d))
    for i in range(d):
        for j in range(d):
            for l in range(d*d):
                if l == (d*i + (i+j)%d):
                    mat[d*i+j][l] = 1
                else:
                    mat[d*i+j][l] = 0
    return mat

def genUt(d):
    mat = np.zeros((d*d, d*d))
    for i in range(d):
        for j in range(d):
            for l in range(d*d):
                if l == (d*((i+j)%d) + j):
                    mat[d*i+j][l] = 1
                else:
                    mat[d*i+j][l] = 0
    return mat

def genV(d, k):
    mat = np.zeros((d*d, d*d))
    for i in range(d):
        for j in range(d):
            for l in range(d*d):
                if l == (d*i + (j+k)%d):
                    mat[d*i+j][l] = 1
                else:
                    mat[d*i+j][l] = 0
    return mat

def genW(d, k):
    mat = np.zeros((d*d, d*d))
    for i in range(d):
        for j in range(d):
            for l in range(d*d):
                if l == (d*((i+k)%d) + j):
                    mat[d*i+j][l] = 1
                else:
                    mat[d*i+j][l] = 0
    return mat

def linear_tran(ct, d, U):
    n = d*d
    for i in range(0, n):
        ct = add(ct, cmult(rot(ct, i), (U[i])))
    return ct

def mat_mult(ct_a, ct_b, d):
    Us = genUs(d)
    Ut = genUt(d)
    ct_a0 = linear_tran(ct_a, d, Us)
    ct_b0 = linear_tran(ct_b, d, Ut)

    ct_ab = mult(ct_a0, ct_b0)

    for i in range(d):
        Vk = genV(d, i)
        Wk = genW(d, i)
        ct_ak = linear_tran(ct_a0, d, Vk)
        ct_bk = linear_tran(ct_b0, d, Wk)
        ct_ab = add(ct_ab, mult(ct_ak, ct_bk))
    
    return ct_ab

ct_ab = mat_mult(a_text, b_text, 3)
print(ct_ab)