from scipy.integrate import quad, dblquad
from commonFuncs import e_repr
from pyprelude.FPToolBox import lmap

# multiplicative integral
def simp_region_integral(dist, func, region):
    # is_prob_dist(dist)
    #if region.type != 'simp':
    #    raise 'it is not a simple region'
    
    if region.dim == 2:
        # two dimensional region
        int_f = dblquad

        def f_mult(p, q):
            return func(p, q) * dist(p, q)

        f_mult = flip(f_mult)
    else:
        # one dimensional region
        int_f = quad
        
        def f_mult(p):
            return func(p) * dist(p)

    return int_f(f_mult, *region.repr)[0]


# calculate all formas of: intgral_{region} (qp)^alpha * dist(p,q)
# where when alpha=(0,0) (2-dim) or 0 (1-dim), it is the measure
def poly_region_int(alpha, dist, region):
    res = -1
    # alpha = tuple(alpha)
    if alpha in region.poly_int.keys():
        # if measure has been calculated before
        res = region.poly_int[alpha]
    else:
        if region.type == 'comp':
            # if it is a composed region
            res_lst = lmap(
                    lambda r: poly_region_int(alpha, dist, r), region.repr)
            res = region.comp_f(res_lst)
        else:
            # if it is a simple region
            if region.dim == 2:
                # two dimensional region
                def func(p, q):
                    return (q ** alpha[0]) * (p ** alpha[1])
            else:
                # one dimensional region
                def func(p):
                    return p ** alpha

            res = simp_region_integral(dist, func, region)

        region.poly_int[alpha] = res    # update measure to the object
    return res


####-----------------------------------------------####
####    Auxiliary functions:
####-----------------------------------------------####

# integration over area
# flip function parameters
def flip(func):
    return lambda y, x: func(x, y)



