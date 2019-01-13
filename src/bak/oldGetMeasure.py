# get measure
# the measure function m for set R
# help function for calc measure, m_Map is for tail recursion
def m_r_help(e, f, region, i, j, p1, p2, q1, q2, m_Map):
    if (e, f, region, i, j) in m_Map.keys():
        return m_Map

    if region == 'L' and i == 1 and j == 1:
        m_Map[e, f, region, i, j] = 1/2 * posi(p1 - p2) * posi(q1 - q2)
        return m_Map

    if region == 'L' and i == 1 and j == 2:
        m_Map[e, f, region, i, j] = 1/2 * posi(p2 - p1) * posi(q2 - q1)
        return m_Map

    if region == 'L' and i == 2 and j == 1:
        m0 = m_r_help(e, f, 'L', 1, 2, p1, p2, q1, q2, m_Map)
        m_Map[e, f, region, i, j] = m0[e, f, 'L', 1, 2]
        return m_Map

    if region == 'L' and i == 2 and j == 2:
        m0 = m_r_help(e, f, 'L', 1, 1, p1, p2, q1, q2, m_Map)
        m_Map[e, f, region, i, j] = m0[e, f, 'L', 1, 1]
        return m_Map

    if region == 'K' and i == 1 and j == 1:
        m_Map[e, f, region, i, j] = p1 * q1
        return m_Map

    if region == 'K' and i == 1 and j == 2:
        m_Map[e, f, region, i, j] = p2 * (1 - q1)
        return m_Map

    if region == 'K' and i == 2 and j == 1:
        m_Map[e, f, region, i, j] = (1 - p1) * q2
        return m_Map

    if region == 'K' and i == 2 and j == 2:
        m_Map[e, f, region, i, j] = (1 - p2) * (1 - q2) 
        return m_Map

    if region == 'R':
        m0 = m_r_help(e, f, 'K', i, j, p1, p2, q1, q2, m_Map)
        m1 = m_r_help(e, f, 'L', i, j, p1, p2, q1, q2, m_Map)
        m_Map[e, f, 'R', i, j] = m0[e, f, 'K', i, j] - m1[e, f, 'L', i, j]
        return m_Map 


# function to calc measures for set S
def m_s_help(e, f, i, m_Map):
    if (e, f, i) in m_Map.keys():
        return m_Map[e, f, i]

    m_Map[e, f, i] = 0.5
    return m_Map
