import numpy as np
from collections import OrderedDict
from commonFuncs import *

#---------------------------------------------------------#
#  main functions to calc E[D]
#---------------------------------------------------------#

# calc expectation from an infomation augmented graph g
# that is run basicInfo on g first
def E_D(g):
    return np.dot(np.dot(g.vec_px, mat_B(g)), g.vec_py)


# calc matrix B
def mat_B(g):

    def f_diag(g, e, f):
        return np.dot(np.dot(np.ones(4), mat_G(g,e,f) * mat_F(g,e,f)), np.ones(3))

    def f_other(g, e, f):
        return np.dot(np.dot(np.ones(4), mat_C(g,e,f) * mat_E(g,e,f)), np.ones(3))

    return edge_mat(g, f_diag, f_other)


## calc vec_ms
#def vec_ms(g, e, f):
#    tmp_ms = OrderedDict([])
#    for i in g.set2:
#        for j in g.set2:
#            r = g.entries[e, f].info.regions['S', i, j]
#            tmp_ms[i, j] = r.measure()
#
#    return od2arr(tmp_ms)

# calc mat_G
def mat_G(g, e, f):
    G = OrderedDict([])
    info = g.entries[e, f].info
    for i in g.set2:
        for j in g.set2:
            entry1 = info.gamma[i,j]
            entry2 = (-1) ** (i + j) * g.d[e]
            entry3 = (-1) ** (i + j - 1) * g.d[f]
            G[i, j] = [entry1, entry2, entry3]

    return od2arr(G)

# calc mat_F
def mat_F(g, e, f):
    F = OrderedDict([])
    info = g.entries[e, f].info
    for i in g.set2:
        for j in g.set2:
            r = info.regions['S', i, j]
            F[i, j] = [r.polyint((0, 0)), r.polyint((0, 1)), r.polyint((1, 0))]

    return od2arr(F)

## calc vec_m
#def vec_m(g,e,f):
#    tmp_m = OrderedDict([])
#
#    for i in g.set2:
#        for j in g.set2:
#            r = g.entries[e, f].info.regions['R', i, j]
#            tmp_m[i, j] = r.measure()
#
#    return od2arr(tmp_m)

# calc mat_C
def mat_C(g, e, f):
    tmp_C = OrderedDict([])
    for i in g.set2:
        for j in g.set2:
            entry1 = g.entries[e,f].info.cs[i,j]
            entry2 = (-1)**(i-1) * g.d[e]
            entry3 = (-1)**(j-1) * g.d[f]
            tmp_C[i,j] = [entry1, entry2, entry3]

    return od2arr(tmp_C)

# calc mat_E
def mat_E(g,e,f):
    E = OrderedDict([])
    info = g.entries[e, f].info

    for i in g.set2:
        for j in g.set2:
            r = info.regions['R', i, j]
            tmpRow = OrderedDict([])
            for alpha in A_curl(2, 1):
                tmpRow[alpha] = r.polyint(alpha)

            E[i, j] = tmpRow

    return od2arr(E)





#    W = mat_W(g,e,f)
#    for alpha in A_curl(2, 1):
#        E[tuple(alpha)] = np.dot(mat_V(g,e,f,alpha) * W, np.ones(2))
#
#    return np.transpose(od2arr(E))


# calc mat_V
#def mat_V(g, e, f, alpha):
#    V = OrderedDict([])
#    info = g.entries[e,f].info
#
#    for i in g.set2:
#        for j in g.set2:
#            V[i,j] = [info.exp(('K', i, j), alpha), info.exp(('L', i, j), alpha)] 
#
#    return od2arr(V)
#
#    
## calc mat_W
#def mat_W(g, e, f):
#    info = g.entries[e,f].info
#    W = OrderedDict([])
#    for i in g.set2:
#        for j in g.set2:
#            denom = info.m['R', i, j]
#            if denom == 0:
#                entry1 = 0
#                entry2 = 0
#            else:
#                entry1 = info.m['K',i,j]/denom
#                entry2 = - info.m['L',i,j]/denom
#            W[i,j] = [entry1, entry2]
#    
#    return od2arr(W)

#---------------------------------------------------------#
#  testing code
#---------------------------------------------------------#

if __name__ == '__main__':

    # od0 = OrderedDict([('c0',0),('c1',1)])
    # od1 = OrderedDict([('c0',2),('c1',3)])
    ls0 = [1, 2, 3]
    ls1 = [4, 5 ,6]
    odd = OrderedDict([('r0',ls0),('r1',ls1)])
    print(ls0, ls1, odd)
    print(od2arrh(ls0))
    print(od2arrh(ls1))
    print(od2arrh(odd))
    
    print(B_curl(3, 2))
    print(A_curl(3, 2)) 
