import networkx as nx
from commonFuncs import *
from collections import OrderedDict
from entryInfo import *
import copy

# if distIndept is true, distInfo contains f_p and f_q
# if distIndept is false, distInfo contains f_pq and f_q|p
def basicInfo(g, distInfo, distIndept=True):


    def uni(p):
        if p >= 0 and p <= 1:
            return 1
        else:
            return 0

    if distInfo == 'uniform':
        g.dist = lambda p, q: uni(p) * uni(q)
        g.distq_p = lambda p, q: uni(q)
    elif distIndept is True:
        distp = distInfo[0]
        distq = distInfo[1]
        # check if it is a distribution
        is_prob_dist(distp, dim=1)
        is_prob_dist(distq, dim=1)
        # create joint pdf
        g.dist = lambda p, q: distp(p) * distq(q)
        g.distq_p = lambda p, q: distq(q)
    else:
        g.dist = distInfo[0]
        # check if it is a distribution
        is_prob_dist(g.dist)
        g.distq_p = distInfo[1]
        
    # basic set
    set2 = [1, 2]
    
    # get ledge length
    d = {}
    for e in g.edges():
        d[e] = g.edges[e]['d']

    # the map from e v id to nat
    # e_idx, v_idx = get_EV_maps(g)

    # shortest path matrix
    map_D = get_map_D(g)

    # update g info
    g.set2 = set2
    g.d = d
    g.map_D = map_D

    # get all entry infos
    # the max possible distance in a network
    max_d = 0

    entryMap = {}
    # diagnal info
    for e in g.edges():
        e = e_repr(e)
        tmpInfo = entry(g, e, e)
        info = tmpInfo.info
        entryMap[e,e] = tmpInfo
        if info.ub > max_d:
            max_d = info.ub


    # off diagnal info
    for e in g.edges():
        for f in g.edges():
            e = e_repr(e)
            f = e_repr(f)
            tmpInfo = entry(g, e, f)
            info = tmpInfo.info
            entryMap[e, f] = tmpInfo
            if info.ub > max_d:
                max_d = info.ub
                    # print(info.part)
                    # print(e, f)
                    # print(info.b)

    g.max_d = max_d
    g.entries = entryMap
    g.vec_px = vec_p(g, 'x')
    g.vec_py = vec_p(g, 'y')

    return g


#---------------------------------------------------------#
#  supplementary functions for getting basic info
#---------------------------------------------------------#

# get edges and nodes maps
#def get_EV_maps(g):
#    e_idx = {}
#    for idx, e in enumerate(g.edges()):
#        e_idx[e] = idx
#
#    v_idx = {}
#    for idx, v in enumerate(g.nodes()):
#        v_idx[v] = idx
#
#    return (e_idx, v_idx)

# get shortest path matrix D
def get_map_D(g):
    return dict(nx.all_pairs_dijkstra_path_length(g, weight = 'd'))


# get probability vector
def vec_p(g, key):
    p = OrderedDict([])
    for e in g.edges():
        p[e] = g.edges[e][key]

    return od2arr(p)


