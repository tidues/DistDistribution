from commonFuncs import *
import measureG as mg
import measureS as ms
import region as rg
import pdfRegion as pr

class entry:
    def __init__(self, g, e, f):
        self.e = e
        self.f = f
        if e_repr(e) == e_repr(f):
            self.info = diagInfo(g, e)
        else:
            self.info = offDiagInfo(g, e, f)

class diagInfo:
    def __init__(self, g, e):
        # basic info
        self.part = 'S'
        self.d_min = g.map_D[e[0]][e[1]]
        self.d = g.d[e]
        self.d_avg = (self.d + self.d_min)/2
        self.r = self.d_avg / self.d
        self.ub = self.d_avg
        self.keys = []

        set2 = [1, 2]

        for i in set2:
            for j in set2:
                self.keys.append((self.part, i, j))
        
        def l1(p):
            return p - self.r

        def l1_p(p):
            return max(0, l1(p))

        def l2(p):
            return p + self.r

        def l2_p(p):
            return min(1, l2(p))

        self.gamma = {}
        for i in g.set2:
            for j in g.set2:
                self.gamma[i, j] = (i - 1) * (self.d_min + self.d)


        # domain of x in each region
        self.domx = {}
        self.domx['S', 1, 1] = (0, self.ub)
        self.domx['S', 1, 2] = (0, self.ub)
        self.domx['S', 2, 1] = (self.d_min, self.ub)
        self.domx['S', 2, 2] = (self.d_min, self.ub)

        # adjusted x
        self.xn = {}

        for key in self.keys:
            self.xn[key] = zeta(*self.domx[key])

        # S_{ij}^x functions info
        self.regions_x = {}

        def S11(x):
            key = ('S', 1, 1)
            xn = self.xn[key](x)

            def l11(p):
                return p - xn/self.d

            def l11_p(p):
                return max(0, l11(p))

            return rg.Region(
                        ('S', 1, 1),
                        (0, 1, l11_p, lambda p: p), g.dist, g.distq_p)

        def S12(x):
            key = ('S', 1, 2)
            xn = self.xn[key](x)

            def l12(p):
                return p + xn/self.d

            def l12_p(p):
                return min(1, l12(p))

            return rg.Region(
                        ('S', 1, 2),
                        (0, 1, lambda p: p, l12_p), g.dist, g.distq_p)

        def S21(x):
            key = ('S', 2, 1)
            xn = self.xn[key](x)

            def l21(p):
                return p - 1 + (xn - self.d_min)/self.d

            return rg.Region(
                        ('S', 2, 1),
                        (1-(xn - self.d_min)/self.d, 1, lambda p: 0, l21), g.dist, g.distq_p)

        def S22(x):
            key = ('S', 2, 2)
            param = self.domx[key]
            xn = self.xn[key](x)

            def l22(p):
                return p + 1 - (xn - self.d_min)/self.d

            return rg.Region(
                        ('S', 2, 2),
                        (0, (xn - self.d_min)/self.d, l22, lambda p: 1), g.dist, g.distq_p)

        self.regions_x['S', 1, 1] = S11
        self.regions_x['S', 1, 2] = S12
        self.regions_x['S', 2, 1] = S21
        self.regions_x['S', 2, 2] = S22

        self.regions = {}

        for key in self.keys:
            # regions are spectial case of regions_x where x is at high level
            self.regions[key] = self.regions_x[key](self.domx[key][1])

        # info for calculate pdf
        self.pdf_regions = {}
        for key in self.keys:
            self.pdf_regions[key] = pr.pdfRegion(key, g.dist, g.distq_p,
                                de=self.d,
                                lb=self.domx[key][0],
                                ub=self.domx[key][1],
                                dmin=self.d_min)

class offDiagInfo:
    def __init__(self, g, e, f):
        # basic info
        self.part = 'R'
        self.de = g.d[e]
        self.df = g.d[f]
        self.A = get_A(g.set2, g.map_D, e, f)
        self.p1 = get_p1(self.A, self.de)
        self.p2 = get_p2(self.A, self.de)
        self.q1 = get_q1(self.A, self.df)
        self.q2 = get_q2(self.A, self.df)
        self.l3p = get_l3p(self.A, self.de, self.df)
        self.l4p = get_l4p(self.A, self.de, self.df)
        self.b = b_val(self.A, self.de, self.df)
        self.ub = self.b
        self.cs = {}
        self.alphas = {}

        set2 = [1, 2]

        self.keys = []
        for i in set2:
            for j in set2:
                self.keys.append((self.part, i, j))
        
        for i in g.set2:
            for j in g.set2:
                self.cs[i, j] = c_val(self.A[i, j], self.de, self.df, i, j)
                self.alphas[i, j] = alpha_val(self.A, self.de, self.df, i, j)

        # region info
        p1, p2, q1, q2, l3p, l4p = self.p1, self.p2, self.q1, self.q2, self.l3p, self.l4p

        def l11(p):
            return min(q1, l3p(p))

        def l12(p):
            return max(q1, l4p(p))

        def l21(p):
            return min(q2, l4p(p))

        def l22(p):
            return max(q2, l3p(p))

        self.domx = {}
        for i in g.set2:
            for j in g.set2:
                self.domx['R', i, j] = (self.A[i, j], self.ub)

        # adjusted x
        self.xn = {}

        for key in self.keys:
            self.xn[key] = zeta(*self.domx[key])

        # R_{ij}^x functions info
        self.regions_x = {}

        def R11(x):
            key = ('R', 1, 1)
            xn = self.xn[key](x)

            l11p = (xn - self.A[1, 1])/self.de
            p11 = min(self.p1, l11p)

            def l11(p):
                return (-self.de * p + xn - self.A[1, 1]) / self.df

            def q11(p):
                return min(self.q1, l11(p))

            return rg.Region(
                ('R', 1, 1), (0, p11, lambda p: 0, q11), g.dist, g.distq_p)

        def R12(x):
            key = ('R', 1, 2)
            xn = self.xn[key](x)

            l12p = (xn - self.A[1, 2])/self.de
            p12 = min(self.p2, l12p)

            def l12(p):
                return (self.df + self.de * p - xn + self.A[1, 2])/self.df

            def q12(p):
                return max(self.q1, l12(p))

            return rg.Region(
                ('R', 1, 2), (0, p12, q12, lambda p: 1), g.dist, g.distq_p)

        def R21(x):
            key = ('R', 2, 1)
            xn = self.xn[key](x)

            l21p = (self.de - xn + self.A[2, 1])/self.de
            p21 = max(self.p1, l21p)

            def l21(p):
                return (-self.de + self.de * p + xn - self.A[2, 1])/self.df

            def q21(p):
                return min(self.q2, l21(p))

            return rg.Region(
                ('R', 2, 1), (p21, 1, lambda p: 0, q21), g.dist, g.distq_p)

        def R22(x):
            key = ('R', 2, 2)
            xn = self.xn[key](x)

            l22p = (self.de - xn + self.A[2, 2])/self.de
            p22 = max(self.p2, l22p)

            def l22(p):
                return (self.de + self.df
                        - self.de * p - xn + self.A[2, 2])/self.df

            def q22(p):
                return max(self.q2, l22(p))

            return rg.Region(
                ('R', 2, 2), (p22, 1, q22, lambda p: 1), g.dist, g.distq_p)

        self.regions_x['R', 1, 1] = R11
        self.regions_x['R', 1, 2] = R12
        self.regions_x['R', 2, 1] = R21
        self.regions_x['R', 2, 2] = R22

        self.regions = {}

        for key in self.keys:
            # regions are spectial case of regions_x where x is at high level
            self.regions[key] = self.regions_x[key](self.domx[key][1])

        # info for calculate pdf
        self.pdf_regions = {}
        for key in self.keys:
            self.pdf_regions[key] = pr.pdfRegion(key, g.dist, g.distq_p,
                                de=self.de,
                                df=self.df,
                                lb=self.domx[key][0],
                                ub=self.domx[key][1],
                                p1=self.p1,
                                p2=self.p2,
                                A=self.A)

# ---------------------------------------------------------#
#  supplementary functions for getting basic info
# ---------------------------------------------------------#

# get map_A
def get_A(set2, map_D, e, f):
    A = {}
    for i in set2:
        for j in set2:
            A[i,j] = map_D[e[i-1]][f[j-1]]
    return A

# get p1 p2 q1 q2 l3(p) l4(p)
def get_p1(A, de):
    return (-A[1,1] + A[2,1] + de)/(2 * de)

def get_p2(A, de):
    return (-A[1,2] + A[2,2] + de)/(2 * de)

def get_q1(A, df):
    return (-A[1,1] + A[1,2] + df)/(2 * df)

def get_q2(A, df):
    return (-A[2,1] + A[2,2] + df)/(2 * df)

def get_l3p(A, de, df):
    return lambda p: (-A[1,1] + A[2,2] + de + df)/(2*df) - de/df * p

def get_l4p(A, de, df):
    return lambda p: (A[1,2] - A[2,1] - de + df)/(2*df) + de/df * p


# get value c
def c_val(A_ij, d_e, d_f, i, j):
    return A_ij + posi((-1) ** i) * d_e + posi((-1) ** j) * d_f

# get alpha i j
def alpha_val(A, d_e, d_f, i, j):
    if (i, j) == (1, 1):
        res = (A[1,1] + A[1,2] + d_f)/2
    elif (i, j) == (1, 2):
        res = (A[1,2] + A[2,2] + d_e)/2
    elif (i, j) == (2, 1):
        res = (A[1,1] + A[2,1] + d_e)/2
    elif (i, j) == (2, 2):
        res = (A[2,1] + A[2,2] + d_f)/2
    else:
        res = 'error: incorret indices'
    return res

# get b e f
def b_val(A, d_e, d_f):
    return (d_e + d_f + min(A[1,1] + A[2,2], A[1,2] + A[2,1]))/2


