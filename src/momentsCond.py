import numpy as np
from collections import OrderedDict
from commonFuncs import *

#---------------------------------------------------------#
#  main functions to calc arbitrary E[D^k]
#---------------------------------------------------------#

# calc expectation from an infomation augmented graph g
# that is run basicInfo on g first
def momentsCond(k, g, e, p):
    return np.dot(vec_B(k, g, e, p), g.vec_py)

# calc matrix B^k
def vec_B(k, g, e, p):

    def f_diag(g, f, k, p):
        return np.dot(np.dot(np.ones(4), mat_G(k, g, e, f) * mat_F(k, g, e, f, p)), np.ones(A_curl_num(2, k)))

    def f_other(g, f, k, p):
        return np.dot(np.dot(np.ones(4), mat_C(k, g, e, f) * mat_E(k, g, e, f, p)), np.ones(A_curl_num(2, k)))
    
    return edge_vec(g, e, f_diag, f_other, k, p)

def mat_G(k, g, e, f):
    G = OrderedDict([])
    info = g.entries[e,f].info
    for i in g.set2:
        for j in g.set2:
            row = OrderedDict([])
            for alpha in A_curl(2, k):
                w = (i + j -1) * alpha[0] + (i + j) * alpha[1]
                row[alpha] = ((-1) ** w) * mChoose(k, alpha) * (g.d[e] ** (alpha[0] + alpha[1])) * (info.gamma[i, j] ** (k - alpha[0] - alpha[1]))

            G[i,j] = row

    return od2arr(G)


def mat_F(k, g, e, f, p):
    F = OrderedDict([])
    info = g.entries[e, f].info
    for i in g.set2:
        for j in g.set2:
            r = info.regions['S', i, j]
            row = OrderedDict([])
            for alpha in A_curl(2, k):
                #print(('S', i, j), ':')
                #print(r.regionCond(p))
                row[alpha] = r.polyintCond(alpha, p)

            F[i, j] = row

    return od2arr(F)

def mat_C(k, g, e, f):
    C = OrderedDict([])
    info = g.entries[e,f].info
    for i in g.set2:
        for j in g.set2:
            row = OrderedDict([])
            c = info.cs[i,j]
            for alpha in A_curl(2, k):
                power = (i - 1) * alpha[1] + (j - 1) * alpha[0]
                row[alpha] = ((-1) ** power) * mChoose(k, alpha) * (g.d[f] ** alpha[0]) * (g.d[e] ** alpha[1]) * (c ** (k - alpha[0] - alpha[1]))
            C[i,j] = row

    return od2arr(C)

def mat_E(k, g, e, f, p):
    E = OrderedDict([])
    info = g.entries[e, f].info

    for i in g.set2:
        for j in g.set2:
            r = info.regions['R', i, j]
            row = OrderedDict([])
            for alpha in A_curl(2, k):
                row[alpha] = r.polyintCond(alpha, p)

            E[i, j] = row

    return od2arr(E)
