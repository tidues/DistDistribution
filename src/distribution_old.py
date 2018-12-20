import numpy as np
from collections import OrderedDict
from commonFuncs import *
from pyprelude.FPToolBox import *

# pdf
def f_D(g):
    def f_D0(x):
        return np.dot(np.dot(g.vec_px, mat_F(g)(x)), g.vec_py)
    return f_D0

def mat_F(g):

    def mat_F0(x):
        def f_diag(g, e, f, x):
            info = g.entries[e,f].info
            ub = info.d_avg
            d_min = info.d_min
            return 2/(g.d[e] ** 2) * (eta(0, ub)(x) * (g.d[e] - x) + eta(d_min, ub)(x) * (x - info.d_min))
        
        def f_other(g, e, f, x):
            info = g.entries[e,f].info
            A = info.A
            b = info.b
            alphas = info.alphas
            term0 = 1/(g.d[e] * g.d[f])
            term1 = 0
            for i in g.set2:
                for j in g.set2:
                    term1 += (eta(A[i,j], b)(x) * (x - A[i,j]) 
                            - 2 * eta(alphas[i,j], b)(x) * (x - alphas[i,j]))
            return term0 * term1

        return edge_mat(g, f_diag, f_other, x)

    return mat_F0

        
# CDF
def F_D(g):
    def F_D0(x):
        return np.dot(np.dot(g.vec_px, mat_CDF(g)(x)), g.vec_py)
    return F_D0


def mat_CDF(g):
    def mat_F0(x):
        def f_diag(g, e, f, x):
            info = g.entries[e,f].info
            ub = info.d_avg
            d_min = info.d_min
            f1 = delta(0, ub, lambda x: 2 * g.d[e] * x - x ** 2)
            f2 = delta(d_min, ub, lambda x: (x - d_min) ** 2)
            return 1/(g.d[e] ** 2) * (f1(x) + f2(x))
        
        def f_other(g, e, f, x):
            info = g.entries[e,f].info
            A = info.A
            b = info.b
            alphas = info.alphas
            term0 = 1/(2 * g.d[e] * g.d[f])
            term1 = 0
            for i in g.set2:
                for j in g.set2:
                    f1 = delta(A[i,j], b, lambda x: (x - A[i,j]) ** 2)
                    f2 = delta(alphas[i,j], b, lambda x: (x - alphas[i,j]) ** 2)
                    term1 += f1(x) - 2 * f2(x)
                            
            return term0 * term1

        return edge_mat(g, f_diag, f_other, x)
    return mat_F0

