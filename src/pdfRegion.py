from commonFuncs import *
from scipy.integrate import quad

class pdfRegion:
    def __init__(self, key, dist, distq_p, **params):
        self.key = key
        self.dist = dist
        self.distq_p = distq_p
        if key[0] == 'S':
            de = params['de']
            lb = params['lb']
            ub = params['ub']
            dmin = params['dmin']
            df = de
        elif key[0] == 'R':
            de = params['de']
            df = params['df']
            lb = params['lb']
            ub = params['ub']
            p1 = params['p1']
            p2 = params['p2']
            A = params['A']
        else:
            raise 'pdfRegion_init_: wrong key 01'
            
        zab = zeta(lb, ub)
        eab = eta(lb, ub)
        self.dxq = lambda x: eab(x)/df

        if key == ('S', 1, 1):
            u = lambda x: zab(x)/de
            self.L = lambda x: (u(x), 1)
            self.q = lambda x, p: p - u(x)
        elif key == ('S', 1, 2):
            u = lambda x: zab(x)/de
            self.L = lambda x: (0, 1 - u(x))
            self.q = lambda x, p: p + u(x)
        elif key == ('S', 2, 1):
            u = lambda x: (zab(x) - dmin)/de
            self.L = lambda x: (1 - u(x), 1)
            self.q = lambda x, p: p - 1 + u(x)
        elif key == ('S', 2, 2):
            u = lambda x: (zab(x) - dmin)/de
            self.L = lambda x: (0, u(x))
            self.q = lambda x, p: p + 1 - u(x)
        elif key == ('R', 1, 1):
            self.L = lambda x: (
                        (2 * zab(x) - (A[1, 1] + A[1, 2] + df))/(2 * de),
                        min(p1, (zab(x) - A[1, 1])/de)
                        )
            self.q = lambda x, p: (-de*p - A[1, 1] + zab(x))/df
        elif key == ('R', 1, 2):
            self.L = lambda x: (
                        (2 * zab(x) - (A[1, 1] + A[1, 2] + df))/(2 * de),
                        min(p2, (zab(x) - A[1, 2])/de)
                        )
            self.q = lambda x, p: (de*p + df + A[1, 2] - zab(x))/df
        elif key == ('R', 2, 1):
            self.L = lambda x: (
                        max(p1, (-zab(x) + A[2, 1] + de)/de),
                        (-2*zab(x) + A[2, 1] + A[2, 2] + 2*de + df)/(2 * de)
                        )
            self.q = lambda x, p: (de*(1-p) - A[2, 1] + zab(x))/df
        elif key == ('R', 2, 2):
            self.L = lambda x: (
                        max(p2, (-zab(x) + A[2, 2] + de)/de),
                        (-2*zab(x) + A[2, 1] + A[2, 2] + 2*de + df)/(2 * de)
                        )
            self.q = lambda x, p: (de*(1-p) + df + A[2, 2] - zab(x))/df
        else:
            raise 'pdfRegion_init_: wrong key 02'

    def pdf(self, x):
        if self.dxq(x) == 0:
            return 0

        def tmpf(p):
            return self.dist(p, self.q(x, p)) * self.dxq(x)

        region = self.L(x)

        return quad(tmpf, *region)[0]

    def pdfCond(self, p):
        def pdf_help(x):
            if self.dxq(x) == 0:
                return 0

            region = self.L(x)
            
            if p >= region[0] and p < region[1]:
                return self.dxq(x) * self.distq_p(p, self.q(x, p))
            else:
                return 0
        return pdf_help







