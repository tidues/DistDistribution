import numpy as np
import pyprelude.FPToolBox as fp
import operator as op
from functools import reduce
from scipy.integrate import quad
from collections import OrderedDict
from scipy.integrate import quad, dblquad

# get max(0,x)
def posi(x):
    if x > 0:
        return x

    return 0

# get unique representation of an edge
def e_repr(e):
    return tuple(sorted(e))

# transform ordereddict to numpy array
def od2arr(od):
    return np.array(od2arrh(od))

# help function of od2arr 
def od2arrh(od):
    # the basic case
    if 'OrderedDict' not in str(type(od)):
        return od

    # recursive call
    return list(map(od2arrh, list(od.values())))

# return ordered set for B_curl(a,b)
def B_curl(n, k):
    tmp = B_curl_h(n, k)
    return list(map(tuple, tmp))

def B_curl_h(n, k):
    if n == 1:
        return [[k]]

    if k == 0:
        return fp.lmap(lambda xs: addhd(0, xs), B_curl_h(n - 1, 0))
    
    # for each possible value, get one
    res = []
    for i in range(k + 1):
        res += fp.lmap(lambda xs: addhd(i, xs), B_curl_h(n - 1, k - i))

    return res

# the number of elements in B_curl
def B_curl_num(n, k):
    if n == 1:
        return 1

    if k == 0:
        return 1
    
    res = 0
    for i in range(k + 1):
        res += B_curl_num(n - 1, k - i)

    return res

# return ordered set for A_curl(a,b)
def A_curl(n, k):
    res = []
    for k0 in range(k + 1):
        res += B_curl(n, k0)

    return res

# return ordered set for A_curl(a,b)
def A_curl_num(n, k):
    res = 0
    for k0 in range(k + 1):
        res += B_curl_num(n, k0)

    return res

# function to fast add head to a list
def addhd(hd, xs):
    xs.insert(0, hd)
    return xs


# value similarity check: values are same with same keys
def similarMaps(map1, map2, epsilon):
    errMsg = 'all correct'
    errMap = {}
    for key in map1.keys():
        if key in map2.keys():
            if abs(map1[key] - map2[key]) > epsilon:
                errMsg = 'something wrong'
                errMap[key] = (map1[key], map2[key])

    return (errMsg, errMap)


# choose function
def choose(n, r):
    r = min(r, n-r)
    numer = reduce(op.mul, range(n, n-r, -1), 1)
    denom = reduce(op.mul, range(1, r+1), 1)
    return numer//denom

# multichoose function
def mChoose(n, alpha):
    x = 1
    for a in alpha:
        x *= choose(n, a)
        n -= a
        if n < 0:
            return 'error: |alpha| > n'
    return x

# the eta function
def eta(a, b):
    def eta1(x):
#        if a is None:
#            a = - float('inf')
#
#        if b is None:
#            b = float('inf')

        if x >= a and x <= b:
            return 1
        else:
            return 0
    return eta1

# the zeta function
def zeta(a, b):
    def zeta1(x):
        if x < a:
            return a
        elif x > b:
            return b
        else:
            return x
    return zeta1

# delta functor
def delta(a, b, f):
    def delta1(x):
        if x < a:
            return f(a)
        elif x > b:
            return f(b)
        else:
            return f(x)
    return delta1

# 1-dim distribution verify
def pdfVeri(func, lb, ub):
    errMsg = ''
    # check if integration over lb ub equal to 1
    total0 = quad(func, lb, ub)
    if contain(total0, 1) is False:
        errMsg = 'not a distribution in lb to ub:\t' + str(total0)
        return errMsg

    # check if intergration over all equal to 1
    # total1 = quad(func, -np.inf, np.inf)
    delta = ub - lb
    total1 = quad(func, lb - delta, ub)
    if contain(total1, 1) is False:
        errMsg = 'not a distribution in whole space:\t' + str(total1)
        return errMsg

    # check if to large difference
    diff = abs(total1[0] - total0[0])
    err = min(total0[1], total1[1])
    if diff > err:
        errMsg = ('may exist nonzero value outside lb to ub:\t' 
                + str(total0) + '\t' + str(total1))
        return errMsg

    return 'all correct, diff:\t' + str(diff)


# check if value is in the range
def contain(vTuple, val):
    res, err = vTuple
    if val < res - err or val > res + err: 
        return False
    else:
        return True

# make edge matrix
def edge_mat(g, f_diag, f_other, *args):
    res = OrderedDict([])
    for e in g.edges():
        tmpRow = OrderedDict([])
        for f in g.edges():
            if e_repr(e) == e_repr(f):
                tmpRow[f] = f_diag(g, e, f, *args)
            else:
                tmpRow[f] = f_other(g, e, f, *args)
        res[e] = tmpRow

    return od2arr(res)

def edge_vec(g, e, f_diag, f_other, *args):
    res = OrderedDict([])
    for f in g.edges():
        if e_repr(e) == e_repr(f):
            res[f] = f_diag(g, f, *args)
        else:
            res[f] = f_other(g, f, *args)

    return od2arr(res)

# test if integration is to 1
def is_prob_dist(func, dim=2):
    if dim ==2:
        res, err = dblquad(func, 0, 1, lambda p: 0, lambda p: 1)
    else:
        res, err = quad(func, 0, 1)

    if 1 < res - err or 1 > res + err:
        raise Exception('error: function provided does not integrate to 1')

    return True
