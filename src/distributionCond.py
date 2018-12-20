import numpy as np
from collections import OrderedDict
from commonFuncs import *
from pyprelude.FPToolBox import *
from scipy.misc import derivative

# CDF
def F_D(g, e, p):
    def F_D0(x):
        return np.dot(vec_CDF(g, e, p)(x), g.vec_py)
    return F_D0

def vec_CDF(g, e, p):
    def vec_F0(x):
        def f_diag(g, f, p, x):
            info = g.entries[e, f].info
            mx = 0
            for i in g.set2:
                for j in g.set2:
                    mx += info.regions_x['S', i, j](x).measureCond(p)

            return mx
        
        def f_other(g, f, p, x):
            info = g.entries[e, f].info
            mx = 0
            for i in g.set2:
                for j in g.set2:
                    mx += info.regions_x['R', i, j](x).measureCond(p)
                            
            return mx

        return edge_vec(g, e, f_diag, f_other, p, x)
    return vec_F0

# pdf by numeric methods
# def f_D1(g, e, p):
#     def f_D0(x):
#         return derivative(F_D(g, e, p), x, dx=1e-6)
#     return f_D0

# pdf
def f_D(g, e, p):
    def f_D0(x):
        return np.dot(vec_PDF(g, e, p)(x), g.vec_py)
    return f_D0

def vec_PDF(g, e, p):
    def vec_f0(x):
        def f_diag(g, f, p, x):
            info = g.entries[e, f].info
            mx = 0
            for i in g.set2:
                for j in g.set2:
                    mx += info.pdf_regions['S', i, j].pdfCond(p)(x)

            return mx
        
        def f_other(g, f, p, x):
            info = g.entries[e, f].info
            mx = 0
            for i in g.set2:
                for j in g.set2:
                    mx += info.pdf_regions['R', i, j].pdfCond(p)(x)
                            
            return mx

        return edge_vec(g, e, f_diag, f_other, p, x)
    return vec_f0
