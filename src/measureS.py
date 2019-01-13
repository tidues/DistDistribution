import numpy as np
from pyprelude.FPToolBox import *
from commonFuncs import *

# function that describe all boundary of regions
def measure_s(re):
    m_Map = {}
    m_Map['T', 1] = 1 / 2
    m_Map['T', 2] = 1 / 2
    m_Map['S', 1, 1] = re - re ** 2 / 2
    m_Map['S', 1, 2] = m_Map['S', 1, 1]
    m_Map['S', 2, 1] = (1 - re) ** 2 / 2
    m_Map['S', 2, 2] = m_Map['S', 2, 1]
    return m_Map

def measure_r(p1, p2, q1, q2):
    m_Map = {}
    set2 = [1, 2]

    m_Map['L', 1, 1] = 1/2 * posi(p1 - p2) * posi(q1 - q2)
    m_Map['L', 1, 2] = 1/2 * posi(p2 - p1) * posi(q2 - q1)
    m_Map['L', 2, 1] = m_Map['L', 1, 2]
    m_Map['L', 2, 2] = m_Map['L', 1, 1]

    m_Map['K', 1, 1] = p1 * q1
    m_Map['K', 1, 2] = p2 * (1 - q1)
    m_Map['K', 2, 1] = (1 - p1) * q2
    m_Map['K', 2, 2] = (1 - p2) * (1 - q2)

    for i in set2:
        for j in set2:
            m_Map['R', i, j] = m_Map['K', i, j] - m_Map['L', i, j]

    return m_Map

# expectation for diagonal region
def exp_s(re):
    m_Map = {}
    I = (0,0)
    P = (0,1)
    Q = (1,0)

    m_Map[('S', 1, 1), I] = 1
    m_Map[('S', 1, 1), P] = (3 - re ** 2) / (6 - 3 * re)
    m_Map[('S', 1, 1), Q] = 1 - m_Map[('S', 1, 1), P]
    #m_Map[('S', 1, 1), PP] = None
    #m_Map[('S', 1, 1), PQ] = None
    #m_Map[('S', 1, 1), QQ] = None

    m_Map[('S', 1, 2), I] = 1
    m_Map[('S', 1, 2), Q] = m_Map[('S', 1, 1), P]
    m_Map[('S', 1, 2), P] = 1 - m_Map[('S', 1, 2), Q]
    #m_Map[('S', 1, 2), PP] = None
    #m_Map[('S', 1, 2), PQ] = None
    #m_Map[('S', 1, 2), QQ] = None

    m_Map[('S', 2, 1), I] = 1
    m_Map[('S', 2, 1), P] = (2 + re) / 3
    m_Map[('S', 2, 1), Q] = (1 - re) / 3
    #m_Map[('S', 2, 1), PP] = None
    #m_Map[('S', 2, 1), PQ] = None
    #m_Map[('S', 2, 1), QQ] = None

    m_Map[('S', 2, 2), I] = 1
    m_Map[('S', 2, 2), P] = (1 - re) / 3
    m_Map[('S', 2, 2), Q] = (2 + re) / 3
    #m_Map[('S', 2, 2), PP] = None
    #m_Map[('S', 2, 2), PQ] = None
    #m_Map[('S', 2, 2), QQ] = None

    return m_Map

# calculate expectation over a region
def exp_r(p1, p2, q1, q2):
    m_Map = {}
    I = (0,0)
    P = (0,1)
    Q = (1,0)
    PP = (0,2)
    QQ = (2,0)
    PQ = (1,1)
    p_vec = np.array([p1, p2])
    q_vec = np.array([q1, q2])
    pp_vec = np.array([p1 ** 2, p1 * p2, p2 ** 2])
    pq_vec = np.array([p1 * q1, p1 * q2, p2 * q1, p2 * q2])
    qq_vec = np.array([q1 ** 2, q1 * q2, q2 ** 2])

    vec12 = np.array([1/3, 2/3])
    vec21 = np.array([2/3, 1/3])
    vec123 = np.array([1/6, 2/6, 3/6])
    vec321 = np.array([3/6, 2/6, 1/6])
    vec5331 = np.array([5/12, 3/12, 3/12, 1/12])
    vec3153 = np.array([3/12, 1/12, 5/12, 3/12])
    vec3513 = np.array([3/12, 5/12, 1/12, 3/12])
    vec1335 = np.array([1/12, 3/12, 3/12, 5/12])

    m_Map[('L', 1, 1), I] = 1
    m_Map[('L', 1, 1), P] = np.dot(vec21, p_vec)
    m_Map[('L', 1, 1), Q] = np.dot(vec21, q_vec)
    m_Map[('L', 1, 1), PP] = np.dot(vec321, pp_vec)
    m_Map[('L', 1, 1), PQ] = np.dot(vec5331, pq_vec)
    m_Map[('L', 1, 1), QQ] = np.dot(vec321, qq_vec)

    m_Map[('L', 1, 2), I] = 1
    m_Map[('L', 1, 2), P] = np.dot(vec12, p_vec)
    m_Map[('L', 1, 2), Q] = np.dot(vec21, q_vec)
    m_Map[('L', 1, 2), PP] = np.dot(vec123, pp_vec)
    m_Map[('L', 1, 2), PQ] = np.dot(vec3153, pq_vec)
    m_Map[('L', 1, 2), QQ] = np.dot(vec321, qq_vec)

    m_Map[('L', 2, 1), I] = 1
    m_Map[('L', 2, 1), P] = np.dot(vec21, p_vec)
    m_Map[('L', 2, 1), Q] = np.dot(vec12, q_vec)
    m_Map[('L', 2, 1), PP] = np.dot(vec321, pp_vec)
    m_Map[('L', 2, 1), PQ] = np.dot(vec3513, pq_vec)
    m_Map[('L', 2, 1), QQ] = np.dot(vec123, qq_vec)

    m_Map[('L', 2, 2), I] = 1
    m_Map[('L', 2, 2), P] = np.dot(vec12, p_vec)
    m_Map[('L', 2, 2), Q] = np.dot(vec12, q_vec)
    m_Map[('L', 2, 2), PP] = np.dot(vec123, pp_vec)
    m_Map[('L', 2, 2), PQ] = np.dot(vec1335, pq_vec)
    m_Map[('L', 2, 2), QQ] = np.dot(vec123, qq_vec)

    m_Map[('K', 1, 1), I] = 1
    m_Map[('K', 1, 1), P] = p1 / 2
    m_Map[('K', 1, 1), Q] =  q1 / 2
    m_Map[('K', 1, 1), PP] = p1 ** 2 / 3
    m_Map[('K', 1, 1), PQ] = p1 * q1 / 4
    m_Map[('K', 1, 1), QQ] = q1 ** 2 / 3

    m_Map[('K', 1, 2), I] = 1
    m_Map[('K', 1, 2), P] = p2 / 2
    m_Map[('K', 1, 2), Q] = (q1 + 1) / 2
    m_Map[('K', 1, 2), PP] = p2 ** 2 / 3
    m_Map[('K', 1, 2), PQ] = p2 * (q1 + 1) / 4
    m_Map[('K', 1, 2), QQ] = (q1 ** 2 + q1 + 1) / 3

    m_Map[('K', 2, 1), I] = 1
    m_Map[('K', 2, 1), P] = (p1 + 1) / 2
    m_Map[('K', 2, 1), Q] =  q2 / 2
    m_Map[('K', 2, 1), PP] = (p1 ** 2 + p1 + 1) / 3
    m_Map[('K', 2, 1), PQ] = (p1 + 1) * q2 / 4
    m_Map[('K', 2, 1), QQ] = q2 ** 2 / 3

    m_Map[('K', 2, 2), I] = 1
    m_Map[('K', 2, 2), P] = (p2 + 1) / 2
    m_Map[('K', 2, 2), Q] = (q2 + 1) / 2
    m_Map[('K', 2, 2), PP] = (p2 ** 2 + p2 + 1) / 3
    m_Map[('K', 2, 2), PQ] = (p2 + 1) * (q2 + 1) / 4
    m_Map[('K', 2, 2), QQ] = (q2 ** 2 + q2 + 1) / 3

    return m_Map

