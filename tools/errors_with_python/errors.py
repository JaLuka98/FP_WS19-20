from sympy import *


def error(f, err_vars=None):
    from sympy import Symbol, latex
    s = 0
    latex_names = dict()

    if err_vars == None:
        err_vars = f.free_symbols

    for v in err_vars:
        err = Symbol('latex_std_' + v.name)
        s += f.diff(v)**2 * err**2
        latex_names[err] = '\\sigma_{' + latex(v) + '}'

    return latex(sqrt(s), symbol_names=latex_names)


Rmax = var('C_p, alpha, T, kappa, V0, Tp')

f = C_p - 9*alpha**2 * kappa * V0* Tp


print(error(f))

#(b2-b1)/(a1-a2)
