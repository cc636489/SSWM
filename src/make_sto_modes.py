

from chaospy import poly
import sympy as sp
from numpy import pi


def make_sto_modes(pc_basis_str, *args):

    dim = pc_basis_str.get("dim", False)
    orthogonal_pol = pc_basis_str.get("basis", False)
    joint_pdf = pc_basis_str.get("joint_pdf", False)
    coefficient = pc_basis_str.get("co_eff", False)

    if all([dim, orthogonal_pol, joint_pdf, coefficient]) is False:
        raise TypeError("obtain NULL value from pc basis structure! Exit program...")

    poly.base.POWER = "**"
    poly.base.SEP = "*"
    str_x = "x = sp.Symbol('x[0]')"
    str_y = "y = sp.Symbol('x[1]')"
    str_t = "t = sp.Symbol('t')"
    exec str_x in globals()
    exec str_y in globals()
    exec str_t in globals()
    for i in range(dim):
        str_dim = "q"+str(i) + "=sp.Symbol('q"+str(i)+"')"
        exec str_dim in globals()

    # create modes one by one.
    string = "[sp.printing.ccode(  sp.integrate( eval(args[arg]) * (expr) * joint_pdf, "
    for j in range(dim):
        string = string + "(q" + str(j) + ",coefficient[2*"+str(j)+"], coefficient[2*"+str(j)+"+1] )"
        if j < dim-1:
            string = string + ", "
        else:
            string = string + ")  )"
    string = string + " for expr in eval(str(orthogonal_pol))]"

    # return a single list or nested list
    if len(args) == 1:
        arg = 0  # "arg" is needed for eval(string)
        return eval(string)
    else:
        return [eval(string) for arg in range(len(args))]
