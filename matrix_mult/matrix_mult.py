import numpy as np
import math
import sys

np.set_printoptions(threshold=sys.maxsize)

def rot(ct, l):
    return np.roll(ct, -l)

def add(a, b):
    return a + b

def mult(a, b):
    return a * b

def cmult(a, b):
    return a * b

# Method 1: multiply vectors with permutation vectors
def genUs(d, k):
    vec = np.zeros(d*d)
    for l in range(d*d):
        if k >= 0:
            if 0 <= l-d*k and l-d*k < d-k:
                vec[l] = 1
        else:
            if -k <= l-(d+k)*d and l-(d+k)*d < d:
                vec[l] = 1
    return vec

def genUt(d, k):
    vec = np.zeros(d*d)
    for i in range(0, d):
        for l in range(d*d):
            if l == k+d*i:
                vec[l] = 1
    return vec 

def genVk(d, k):
    vec = np.zeros(d*d)
    for l in range(d*d):
        if l%d >= 0 and l%d < d-k:
            vec[l] = 1
    return vec

def genVkd(d, k):
    vec = np.zeros(d*d)
    for l in range(d*d):
        if l%d >= d-k and l%d < d:
            vec[l] = 1
    return vec 

def linear_tran(ct, d, U_type, k_m):
    if U_type == 0:
        i = 0
        ct_ls = []
        ct_ls.append(cmult(ct, genUs(d, 0)))
        i += 1
        for k in range(1-d, 0):
            ct_ls.append(add(ct_ls[i-1], cmult(rot(ct, k), genUs(d, k))))
            i += 1
        for k in range(1, d):
            ct_ls.append(add(ct_ls[i-1], cmult(rot(ct, k), genUs(d, k))))
            i += 1
        result = ct_ls[i-1]
    elif U_type == 1:
        i = 0
        ct_ls = []
        ct_ls.append(cmult(ct, genUt(d, 0)))
        i += 1
        for k in range(1, d):
            ct_ls.append(add(ct_ls[i-1], cmult(rot(ct, d*k), genUt(d, k))))
            i += 1
        result = ct_ls[i-1]
    elif U_type == 2:
        result = add(cmult(rot(ct, k_m), genVk(d, k_m)), cmult(rot(ct, k_m-d), genVkd(d, k_m)))
    elif U_type == 3:
        result = rot(ct, d*k_m)
    else:
        print('Wrong U_type!')
    return result


def mat_mult(ct_a, ct_b, d):
    ct_a0 = linear_tran(ct_a, d, 0, 0)
    ct_b0 = linear_tran(ct_b, d, 1, 0)

    i = 0
    ab_ls = []
    ab_ls.append(mult(ct_a0, ct_b0))
    i += 1
    for k in range(1, d):
        ct_ak = linear_tran(ct_a0, d, 2, k)
        ct_bk = linear_tran(ct_b0, d, 3, k)
        ab_ls.append(add(ab_ls[i-1], mult(ct_ak, ct_bk)))
        i += 1
    result = ab_ls[i-1]
    return result

def rmat_mult(ct_a, ct_b, d, l):
    ct_a0 = linear_tran(ct_a, d, 0, 0)
    ct_b0 = linear_tran(ct_b, d, 1, 0)

    i = 0
    ab_ls = []
    ab_ls.append(mult(ct_a0, ct_b0))
    i += 1

    for k in range(1, l):
        ct_ak = linear_tran(ct_a0, d, 2, k)
        ct_bk = linear_tran(ct_b0, d, 3, k)
        ab_ls.append(add(ab_ls[i-1], mult(ct_ak, ct_bk)))
        i += 1

    ab = ab_ls[i-1]
    for k in range(0, math.ceil(math.log(d/l, 2))):
        ab_ls.append(add(ab_ls[i-1], rot(ab, l*d*(2**k))))
        i += 1
    result = ab_ls[i-1]
    return result

a = np.matrix([[0, 1, 2], 
              [3, 4, 5], 
              [6, 7, 8]])
ar = np.matrix([[0, 1, 2]])
arr = np.matrix([[0, 1, 2], 
                [0, 1, 2], 
                [0, 1, 2]])
b = np.matrix([[10, 11, 12], 
              [13, 14, 15], 
              [16, 17, 18]])

print("Square matrix multiplication true result:")
a_b = a.dot(b)
print(a_b)

print('Square matrix multiplication result by algorithm: ')
a_text = a.A1
b_text = b.A1
ct_ab = mat_mult(a_text, b_text, 3)
print(ct_ab)

print('Rectangular matrix multiplication true result: ')
ar_b = ar.dot(b)
print(ar_b)

print('Rectangular matrix multiplication result by algorithm: ')
arr_text = arr.A1
ct_arrb = rmat_mult(arr_text, b_text, 3, 1)
print(ct_arrb)

# Method 2: mutiply transformed matrices
def sigma(A):
    d = A.shape[0]
    B = np.zeros((d, d))
    for i in range(d):
        for j in range(d):
            B[i, j] = A[i, (i+j)%d]
    return B

def tau(A):
    d = A.shape[0]
    B = np.zeros((d, d))
    for i in range(d):
        for j in range(d):
            B[i, j] = A[(i+j)%d, j]
    return B

def phi(A):
    d = A.shape[0]
    B = np.zeros((d, d))
    for i in range(d):
        for j in range(d):
            B[i, j] = A[i, (j+1)%d]
    return B

def psi(A):
    d = A.shape[0]
    B = np.zeros((d, d))
    for i in range(d):
        for j in range(d):
            B[i, j] = A[(i+1)%d, j]
    return B

def mat_mult_t(a, b, d):
    x = []
    y = []
    x.append(sigma(a))
    y.append(tau(b))
    a_b = mult(x[0], y[0])
    for i in range(1, d):
        x.append(phi(x[i-1]))
        y.append(psi(y[i-1]))
        a_b = add(a_b, mult(x[i], y[i]))

    return a_b

print('Square matrix multiplication result for verification: ')
print(mat_mult_t(a, b, 3))

