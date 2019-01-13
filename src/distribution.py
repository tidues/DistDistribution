import numpy as np
from collections import OrderedDict
from commonFuncs import *
from pyprelude.FPToolBox import *
from scipy.misc import derivative

# CDF
def F_D(g):
    def F_D0(x):
        return np.dot(np.dot(g.vec_px, mat_CDF(g)(x)), g.vec_py)
    return F_D0

def mat_CDF(g):
    def mat_F0(x):
        def f_diag(g, e, f, x):
            info = g.entries[e, f].info
            mx = 0
            for i in g.set2:
                for j in g.set2:
                    mx += info.regions_x['S', i, j](x).measure()
            return mx
        
        def f_other(g, e, f, x):
            info = g.entries[e, f].info
            mx = 0
            for i in g.set2:
                for j in g.set2:
                    mx += info.regions_x['R', i, j](x).measure()
                            
            return mx

        return edge_mat(g, f_diag, f_other, x)
    return mat_F0

# pdf by numeric methods
# def f_D(g):
#     def f_D0(x):
#         return derivative(F_D(g), x, dx=1e-6)
#     return f_D0

# pdf
def f_D(g):
    def f_D0(x):
        return np.dot(np.dot(g.vec_px, mat_PDF(g)(x)), g.vec_py)
    return f_D0

def mat_PDF(g):
    def mat_f0(x):
        def f_diag(g, e, f, x):
            info = g.entries[e, f].info
            mx = 0
            for i in g.set2:
                for j in g.set2:
                    mx += info.pdf_regions['S', i, j].pdf(x)

            return mx
        
        def f_other(g, e, f, x):
            info = g.entries[e, f].info
            mx = 0
            for i in g.set2:
                for j in g.set2:
                    #if e == ('1', '2') and f == ('2', '3') and i == 2 and j == 1:
                       # print('vv4:')
                       # print('val:',info.pdf_regions['R', i, j].pdf(x))
                    mx += info.pdf_regions['R', i, j].pdf(x)

            return mx

        return edge_mat(g, f_diag, f_other, x)
    return mat_f0
