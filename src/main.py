import basicInfo as bi
import expectation as ep
import pyprelude.GraphGen as gg
from commonFuncs import *
# from measure import check_measure
import moments as mm
import distribution as dis
import distributionCond as dic
from plot import *
from correctnessCheck import *
import momentsCond as mc
from pyprelude.FPToolBox import concat, foldr
from scipy import integrate

fpath = '../data/'
gname = 'g3'
epsilon = 1e-8
plotOffset = 3
dist = 'uniform'

g = gg.GraphFromFile(fpath, gname).G
# gg = bi.basicInfo(g, 'uniform')

def cdist(p):
    if p >= 0 and p <= 1:
        return 1
    else:
        return 0

def ldist_p(p):
    if p >= 0 and p <= 1:
        return 2 * p
    else:
        return 0

def ldist_n(p):
    if p >= 0 and p <= 1:
        return -2 * p + 2
    else:
        return 0

def ddist1(p, q):
    return cdist(p) * cdist(q)

def ddist2(p, q):
    return cdist(p) * ldist_n(q)

def ddist3(p, q):
    return cdist(p) * ldist_p(q)

def ddist4(p, q):
    return ldist_n(p) * ldist_p(q)

def ddist5(p, q):
    return ldist_p(p) * ldist_p(q)

g = bi.basicInfo(g, [ldist_p, ldist_p])
# region_check(gg, epsilon)
# region_x_check(g, ub=False, part='A')
# region_x_check(gg, ub=False)
# print(measure_sum_check(gg, epsilon))
# print(info_checker(g, region_sameas_region_x_check))
#print(measure_agree_check(g, epsilon))
#for e in gg.edges():
#    for f in gg.edges():
#        print(gg.entries[e,f].info.exps)
#print(exp_agree_check(g, epsilon))

#print("set2:\t", g.set2)
#print("e_idx\t", g.e_idx) #print("v_idx:\t", g.v_idx)
#print("map_D:\t", g.map_D)
#print("map_A:\t", g.map_A)
#print("p1:\t", g.p1)
#print("p2:\t", g.p2)
#print("q1:\t", g.q1)
#print("q2:\t", g.q2)
#print("idx2:\t", g.idx2)
#print("idx2s:\t", g.idx2s)
#print("d:\t", g.d)
#print(check_measure_r(g, g.m_r, epsilon))
#print(check_measure_r(g, g.m_r1, epsilon))
#print(check_measure_s(g, g.m_s, epsilon))
#print(check_measure_s(g, g.m_s1, epsilon))
#print(check_measure(g, g.m_rs, epsilon))

# exp = ep.E_D(g)
exp0 = mm.moments(0, g)
exp1 = mm.moments(1, g)
exp2 = mm.moments(2, g)
exp3 = mm.moments(3, g)
var = exp2 - exp1 ** 2

print(foldr(concat('\t'), ['total', exp0, exp1, exp2, exp3, var], ''))

pLst = [0, 0.2, 0.4, 0.6, 0.85, 0.87, 0.89, 0.9, 0.92, 0.94, 0.96]
# pLst = [0]

#e = list(g.edges())[0]
#for p in pLst:
#    expc0 = mc.momentsCond(0, g, e, p)
#    expc1 = mc.momentsCond(1, g, e, p)
#    expc2 = mc.momentsCond(2, g, e, p)
#    varc = expc2 - expc1 ** 2
#    print(foldr(concat('\t'), [p, expc0, expc1, expc2, varc], ''))
#    print(foldr(concat('\t'), [p, expc0], ''))

#
valLst = [-0.1, -0.01, -0.001, -0.0001, 0, 1, 2, 3, -1, 9, 9.5, 9.9, 10]

pdf_D = dis.f_D(g)
cdf_D = dis.F_D(g)

# mat_D = dis.mat_PDF(g)
# 
# print(mat_D(3))

for x in valLst:
    print(x, '\t', pdf_D(x))

#for p in pLst:
#    pdf_t = dic.f_D(g, e, p)
##    pdf_t1 = dic.f_D1(g, e, p)
##    for x in range(10):
##        diff = abs(pdf_t(x) - pdf_t1(x))
##        if diff > 0.01:
##            print(diff)
#
#    cdf_t = dic.F_D(g, e, p)
##    for v in valLst:
##        if v < 0 and pdf_t(v) > 0:
##            print(v, '\t', pdf_t(v))
##
##    for v in valLst:
##        if v < 0 and cdf_t(v) > 0:
##            print(v, '\t', cdf_t(v))
#    plot1d(pdf_t, 0 - plotOffset, g.max_d + plotOffset, 0.2)
#    plot1d(cdf_t, 0 - plotOffset, g.max_d + plotOffset, 0.2)


#
#print(pdfVeri(pdf_D, 0, 60))


#for v in valLst:
#    print(pdf_D(v))
#
#for v in valLst:
#    print(cdf_D(v))
#
print(g.max_d)
#plot1d(pdf_D, 0 - plotOffset, g.max_d + plotOffset, 0.2)
# plot1d(cdf_D, 0 - plotOffset, g.max_d + plotOffset, 0.2)
