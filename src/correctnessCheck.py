from commonFuncs import *

# check 1). if summation of measures equal to 1
def measure_sum_check(g, epsilon):
    errMap = {}
    errMsg = 'all correct'
    for e in g.edges():
        for f in g.edges():
            e = e_repr(e)
            f = e_repr(f)
            if e <= f:
                info = g.entries[e, f].info
                regions = info.regions
                part = info.part

                # check if sum to 1
                total_m = 0.0
                for r in regions.keys():
                    if r[0] == part:
                        total_m += regions[r].measure()

                if total_m < 1.0 - epsilon or total_m > 1.0 + epsilon:
                    return ('something wrong', e, f, part, total_m)

    return 'all correct'

# check 2). if results from two methods are agree
def measure_agree_check(g, epsilon):
    g1 = copy.deepcopy(g)
    g2 = copy.deepcopy(g)

    dist1 = 'uniform'
    dist2 = lambda x, y: 1

    g1 = basicInfo(g1, dist1)
    g2 = basicInfo(g2, dist2)

    for e in g.edges():
        for f in g.edges():
            e = e_repr(e)
            f = e_repr(f)
            info1 = g1.entries[e,f].info
            info2 = g2.entries[e,f].info
            if e <= f:
                m1 = info1.m
                m2 = info2.m
                res = similarMaps(m1, m2, epsilon)
                if res[0] == 'something wrong':
                    return (res, e, f, m1, m2)

    return 'all correct'

# expecation agree checks
def exp_agree_check(g, epsilon):
    g1 = copy.deepcopy(g)
    g2 = copy.deepcopy(g)

    dist1 = 'uniform'
    dist2 = lambda x, y: 1

    g1 = basicInfo(g1, dist1)
    g2 = basicInfo(g2, dist2)

    for e in g.edges():
        for f in g.edges():
            e = e_repr(e)
            f = e_repr(f)
            if e <= f:
                info1 = g1.entries[e,f].info
                info2 = g2.entries[e,f].info
                for i in g.set2:
                    for j in g.set2:
                        region1 = (info1.part, i, j)
                        region2 = (info2.part, i, j)
                        for alpha in A_curl(2,2):
                            exp1 = info1.exp(region1, alpha)
                            exp2 = info2.exp(region2, alpha)
                            if abs(exp1 - exp2) > epsilon:
                                return ('something wrong', e, f, i, j, alpha, exp1, exp2)

    return 'all correct'


# check whether the total measure relative to p is sum to 1
def region_check(g, epsilon):
    for e in g.edges():
        for f in g.edges():
            e = e_repr(e)
            f = e_repr(f)
            if e <= f:
                info = g.entries[e, f].info
                for i in g.set2:
                    for j in g.set2:
                        regions = info.regions
                        print(region_check_h(regions, 0))
                        print(region_check_h(regions, 0.24))
                        print(region_check_h(regions, 0.73))
                        print(region_check_h(regions, 1))


def region_check_h(regions, p):
    resLst = []
    total = 0
    key = None
    for key in regions:
        part = key[0]
        r = regions[key]
        if p >= r.repr[0] and p < r.repr[1]:
            tmp = r.repr[3](p) - r.repr[2](p)
        else:
            tmp = 0

        resLst.append(tmp)
        total += tmp

    if total != 0 and key is not None:
        return (p, part, total, resLst)

# print the region when at the lb and ub
# ub: True, print the region at ub level
#     False, print the region at lb level
# part: specify printed part: 'A' all; 'S'/'R' corresponding region.
def region_x_check(g, ub=True, part='A'):
    for e in g.edges():
        for f in g.edges():
            e = e_repr(e)
            f = e_repr(f)
            if e <= f:
                info = g.entries[e, f].info
                tmpr = {}
                for i in g.set2:
                    for j in g.set2:
                        if part != 'A' and info.part != part:
                            continue
                        key = (info.part, i, j)
                        if ub is True:
                            val = info.domx[key][1]
                        else:
                            val = info.domx[key][0]
                        tmpr[key] = info.regions_x[key](val)
#                if info.part == 'R':
#                    print(info.p1, info.p2, info.q1, info.q2)
                print(region_check_h(tmpr, 0))
                print(region_check_h(tmpr, 0.24))
                print(region_check_h(tmpr, 0.73))
                print(region_check_h(tmpr, 0.999))

# givien a region, check whether regions is the
# special case where x = ub
def region_sameas_region_x_check(info):
    set2 = [1, 2]
    diffLst = []
    for i in set2:
        for j in set2:
            key = (info.part, i, j)
            r1 = info.regions[key]
            r2 = info.regions_x[key](info.domx[key][1])
            diffLst.append(r1.measure() - r2.measure())
    return diffLst


## check correctness of info by an correct checker
def info_checker(g, func):
    for e in g.edges():
        for f in g.edges():
            e = e_repr(e)
            f = e_repr(f)
            if e <= f:
                info = g.entries[e, f].info
                res = func(info)
                if res is not None:
                    print(res)

# check if two measure functions are agree
def check_measure(g, m, epsilon):
    errorMap = {}
    errMsg = 'all correct'
    for e in g.edges():
        for f in g.edges():
            if e_repr(e) != e_repr(f):
                total_m_s = 0
                for i in g.set2:
                    total_m_s += m[e, f, ('S', i)]
                if total_m_s < 1.0 - epsilon or total_m_s > 1.0 + epsilon:
                    errMsg = 'something wrong'
                    errorMap[e, f, ('S', i)] = total_m_s

                total_m_r = 0
                for i in g.set2:
                    for j in g.set2:
                        total_m_r += m[e, f, ('R', i, j)]
