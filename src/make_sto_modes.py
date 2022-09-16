
from subprocess import check_call
from chaospy import poly
import sympy as sp
from numpy import pi


def make_sto_modes(pc_basis_str, *args):

    dim = pc_basis_str.get("dim", False)
    orthogonal_pol = pc_basis_str.get("basis", False)
    joint_pdf = pc_basis_str.get("joint_pdf", False)
    coefficient = pc_basis_str.get("co_eff", False)

    if dim is False or orthogonal_pol is False or joint_pdf is False or coefficient is False:
        raise TypeError("obtain NULL value from pc basis structure! Exit program...")

    poly.base.POWER = "**"
    poly.base.SEP = "*"
    str_x = "x = sp.Symbol('x[0]')"
    str_y = "y = sp.Symbol('x[1]')"
    str_t = "t = sp.Symbol('t')"
    exec(str_x)
    exec(str_y)
    exec(str_t)
    for i in range(dim):
        str_dim = "q"+str(i) + "=sp.Symbol('q"+str(i)+"')"
        exec(str_dim)

    # create modes one by one.
    if len(args) == 1:
        string = "sp.printing.ccode(  sp.integrate( eval(args[0]) * (expr) * joint_pdf, "
    else:
        string = "sp.printing.ccode(  sp.integrate( eval(args[arg]) * (expr) * joint_pdf, "

    for j in range(dim):
        string = string + "(q" + str(j) + ",coefficient[2*"+str(j)+"], coefficient[2*"+str(j)+"+1] )"
        if j < dim-1:
            string = string + ", "
        else:
            string = string + ")  )"

    # prepare orthogonal_pol valid list for sympy integration.
    orth_list = []
    for orth in orthogonal_pol:
        orth_list.append(eval(str(orth)))

    # construct a return: either a single list or nested list
    ret = []
    if len(args) == 1:
        for expr in orth_list:
            ret.append(eval(string))
    else:
        for arg in range(len(args)):
            tmp = []
            for expr in orth_list:
                tmp.append(eval(string))
            ret.append(tmp)

    # return
    return ret
