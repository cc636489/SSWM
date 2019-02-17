from fenics import *
from makestobasis import MakeStoBasis
from makemesh import MakeMesh
from make_bath import make_bath


def test1_deg0_flat():

    distname = "uniform"
    sto_poly_deg = 0
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bathymetry = {"flat": 20.}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    nmodes = basis_struct.get("nmodes")
    string_Q = "["
    for i in range(nmodes-1):
        string_Q = string_Q + "Q,"
    string_Q = string_Q + "Q]"
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    if nmodes == 1:
        P = FunctionSpace(mesh, "CG", 1)
    else:
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))
    make_bath(bathymetry, basis_struct, P)

    return 1


def test2_deg0_vary():

    distname = "uniform"
    sto_poly_deg = 0
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bathymetry = {"vary": '2*x*q0+3*y*q1'}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    nmodes = basis_struct.get("nmodes")
    string_Q = "["
    for i in range(nmodes-1):
        string_Q = string_Q + "Q,"
    string_Q = string_Q + "Q]"
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    if nmodes == 1:
        P = FunctionSpace(mesh, "CG", 1)
    else:
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))

    make_bath(bathymetry, basis_struct, P)

    return 1


def test3_deg0_class():

    distname = "uniform"
    sto_poly_deg = 0
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bathymetry = {"class": ['type1', 400, 600, '5.0', '-3.0*((x-500.0)/100.0)**4*q0 + 6*((x-500.)/100.)**2*q1 + 2.']}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    nmodes = basis_struct.get("nmodes")
    string_Q = "["
    for i in range(nmodes-1):
        string_Q = string_Q + "Q,"
    string_Q = string_Q + "Q]"
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    if nmodes == 1:
        P = FunctionSpace(mesh, "CG", 1)
    else:
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))
    make_bath(bathymetry, basis_struct, P)


    return 1


def test4_degn_flat():

    distname = "uniform"
    sto_poly_deg = 2
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bathymetry = {"flat": 20.}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    nmodes = basis_struct.get("nmodes")
    string_Q = "["
    for i in range(nmodes-1):
        string_Q = string_Q + "Q,"
    string_Q = string_Q + "Q]"
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    if nmodes == 1:
        P = FunctionSpace(mesh, "CG", 1)
    else:
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))
    make_bath(bathymetry, basis_struct, P)


    return 1


def test5_degn_vary():

    distname = "uniform"
    sto_poly_deg = 2
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bathymetry = {"vary": '2*x*q0+3*y*q1'}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    nmodes = basis_struct.get("nmodes")
    string_Q = "["
    for i in range(nmodes-1):
        string_Q = string_Q + "Q,"
    string_Q = string_Q + "Q]"
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    if nmodes == 1:
        P = FunctionSpace(mesh, "CG", 1)
    else:
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))

    make_bath(bathymetry, basis_struct, P)


    return 1


def test6_degn_class():

    distname = "uniform"
    sto_poly_deg = 2
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bathymetry = {"class": ['type1', 400, 600, '5.0', '-3.0*((x-500.0)/100.0)**4*q0 + 6*((x-500.)/100.)**2*q1 + 2.']}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    nmodes = basis_struct.get("nmodes")
    string_Q = "["
    for i in range(nmodes-1):
        string_Q = string_Q + "Q,"
    string_Q = string_Q + "Q]"
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    if nmodes == 1:
        P = FunctionSpace(mesh, "CG", 1)
    else:
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))
    make_bath(bathymetry, basis_struct, P)


    return 1


def main():

    red = '\033[91m'
    green = '\033[92m'
    blue = '\033[94m'
    bold = '\033[1m'

    flag = test1_deg0_flat()
    if flag:
        print blue+bold+"Test1: "+bold+green+"PASS"
    else:
        print blue+bold+"Test1: "+bold+red+"FAIL"

    flag = test2_deg0_vary()
    if flag:
        print blue+bold+"Test2: "+bold+green+"PASS"
    else:
        print blue+bold+"Test2: "+bold+red+"FAIL"

    flag = test3_deg0_class()
    if flag:
        print blue+bold+"Test3: "+bold+green+"PASS"
    else:
        print blue+bold+"Test3: "+bold+red+"FAIL"

    flag = test4_degn_flat()
    if flag:
        print blue+bold+"Test4: "+bold+green+"PASS"
    else:
        print blue+bold+"Test4: "+bold+red+"FAIL"

    flag = test5_degn_vary()
    if flag:
        print blue+bold+"Test5: "+bold+green+"PASS"
    else:
        print blue+bold+"Test5: "+bold+red+"FAIL"

    flag = test6_degn_class()
    if flag:
        print blue+bold+"Test6: "+bold+green+"PASS"
    else:
        print blue+bold+"Test6: "+bold+red+"FAIL"

    return 0


if __name__ == '__main__':
    main()
