import measureG as mg
from scipy.integrate import quad

class Region:
    def __init__(self, ident, rrepr, dist, distq_p, rtype='simp', comp_f=None, dim=2):
        self.dim = dim

        if dim == 2:
            self.zero = (0, 0)
        else:
            self.zero = 0

        self.id = ident          # the region id
        self.type = rtype        # simp: simple or comp: composed
        self.comp_f = comp_f     # the merge function for composed region
        self.repr = rrepr        # the region representation
        # all the moments info, the a_th entry is 
        # E[(QP)^a | region] * m(region)
        # where when a = (0,0), it is the probability
        self.poly_int = {}
        # the conditional moments info
        # E[(QP)^a | region] = poly_int[alpha]/poly_int[0]
        # self.moments = {}
        self.dist = dist
        self.distq_p = distq_p

    def measure(self):
        return self.polyint(self.zero)

    def measureCond(self, p):
        return self.polyintCond(self.zero, p)

    def polyint(self, alpha):
        return mg.poly_region_int(alpha, self.dist, self)

    def regionCond(self, p):
        if p >= self.repr[0] and p < self.repr[1]:
            return (self.repr[2](p), self.repr[3](p))
        else:
            return None


    def polyintCond(self, alpha, p):
        tmpr = self.regionCond(p)

        if tmpr is None:
            return 0

        const = p ** alpha[1]

        def dist(q):
            return q ** alpha[0] * self.distq_p(p, q)

        return quad(dist, *tmpr)[0] * const

