from fenics import *
from makemesh import MakeMesh
from makestobasis import MakeStoBasis
from makeic import MakeInitialObjectFunc

def test1_deg0_flat_flat():

    distname = "uniform"
    sto_poly_deg = 0
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    initial_u = {"flat": (0., 0.)}
    initial_eta = {"flat": 0.0}
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    V = VectorElement(family='CG', cell=mesh.ufl_cell(), degree=2, dim=2)
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    nmodes = basis_struct.get("nmodes")
    string_V = "["
    string_Q = "["
    for i in range(nmodes-1):
        string_V = string_V + "V,"
        string_Q = string_Q + "Q,"
    string_V = string_V + "V]"
    string_Q = string_Q + "Q]"
    if nmodes == 1:
        W = VectorFunctionSpace(mesh, "CG", 2, dim=2)
        P = FunctionSpace(mesh, "CG", 1)
    else:
        W = FunctionSpace(mesh, MixedElement(eval(string_V)))
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))

    MakeInitialObjectFunc(initial_u, initial_eta, basis_struct, W, P)

    return 1


def test2_deg0_flat_vary():

    distname = "uniform"
    sto_poly_deg = 0
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    initial_u = {"flat": (0., 0.)}
    initial_eta = {"vary": '20.*x*y'}
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    V = VectorElement(family='CG', cell=mesh.ufl_cell(), degree=2, dim=2)
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    nmodes = basis_struct.get("nmodes")
    string_V = "["
    string_Q = "["
    for i in range(nmodes-1):
        string_V = string_V + "V,"
        string_Q = string_Q + "Q,"
    string_V = string_V + "V]"
    string_Q = string_Q + "Q]"
    if nmodes == 1:
        W = VectorFunctionSpace(mesh, "CG", 2, dim=2)
        P = FunctionSpace(mesh, "CG", 1)
    else:
        W = FunctionSpace(mesh, MixedElement(eval(string_V)))
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))

    MakeInitialObjectFunc(initial_u, initial_eta, basis_struct, W, P)

    return 1


def test3_deg0_vary_flat():
    distname = "uniform"
    sto_poly_deg = 0
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    initial_u = {"vary": ('2*sp.sin(x+y)', '3*x**4')}
    initial_eta = {"flat": 0.0}
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    V = VectorElement(family='CG', cell=mesh.ufl_cell(), degree=2, dim=2)
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    nmodes = basis_struct.get("nmodes")
    string_V = "["
    string_Q = "["
    for i in range(nmodes-1):
        string_V = string_V + "V,"
        string_Q = string_Q + "Q,"
    string_V = string_V + "V]"
    string_Q = string_Q + "Q]"
    if nmodes == 1:
        W = VectorFunctionSpace(mesh, "CG", 2, dim=2)
        P = FunctionSpace(mesh, "CG", 1)
    else:
        W = FunctionSpace(mesh, MixedElement(eval(string_V)))
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))

    MakeInitialObjectFunc(initial_u, initial_eta, basis_struct, W, P)

    return 1


def test4_deg0_vary_flat():
    distname = "uniform"
    sto_poly_deg = 0
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    initial_u = {"vary": ('2*sp.sin(x+y)', '3*x**4')}
    initial_eta = {"vary": '20.*x*y'}
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    V = VectorElement(family='CG', cell=mesh.ufl_cell(), degree=2, dim=2)
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    nmodes = basis_struct.get("nmodes")
    string_V = "["
    string_Q = "["
    for i in range(nmodes-1):
        string_V = string_V + "V,"
        string_Q = string_Q + "Q,"
    string_V = string_V + "V]"
    string_Q = string_Q + "Q]"
    if nmodes == 1:
        W = VectorFunctionSpace(mesh, "CG", 2, dim=2)
        P = FunctionSpace(mesh, "CG", 1)
    else:
        W = FunctionSpace(mesh, MixedElement(eval(string_V)))
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))

    MakeInitialObjectFunc(initial_u, initial_eta, basis_struct, W, P)

    return 1


def test5_degn_flat_flat():
    distname = "uniform"
    sto_poly_deg = 2
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    initial_u = {"flat": (0., 0.)}
    initial_eta = {"flat": 0.0}
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    V = VectorElement(family='CG', cell=mesh.ufl_cell(), degree=2, dim=2)
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    nmodes = basis_struct.get("nmodes")
    string_V = "["
    string_Q = "["
    for i in range(nmodes-1):
        string_V = string_V + "V,"
        string_Q = string_Q + "Q,"
    string_V = string_V + "V]"
    string_Q = string_Q + "Q]"
    if nmodes == 1:
        W = VectorFunctionSpace(mesh, "CG", 2, dim=2)
        P = FunctionSpace(mesh, "CG", 1)
    else:
        W = FunctionSpace(mesh, MixedElement(eval(string_V)))
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))

    MakeInitialObjectFunc(initial_u, initial_eta, basis_struct, W, P)

    return 1


def test6_degn_flat_vary():
    distname = "uniform"
    sto_poly_deg = 2
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    initial_u = {"flat": (0., 0.)}
    initial_eta = {"vary": '20.*x*y'}
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    V = VectorElement(family='CG', cell=mesh.ufl_cell(), degree=2, dim=2)
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    nmodes = basis_struct.get("nmodes")
    string_V = "["
    string_Q = "["
    for i in range(nmodes-1):
        string_V = string_V + "V,"
        string_Q = string_Q + "Q,"
    string_V = string_V + "V]"
    string_Q = string_Q + "Q]"
    if nmodes == 1:
        W = VectorFunctionSpace(mesh, "CG", 2, dim=2)
        P = FunctionSpace(mesh, "CG", 1)
    else:
        W = FunctionSpace(mesh, MixedElement(eval(string_V)))
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))

    MakeInitialObjectFunc(initial_u, initial_eta, basis_struct, W, P)

    return 1


def test7_degn_vary_flat():
    distname = "uniform"
    sto_poly_deg = 2
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    initial_u = {"vary": ('2*sp.sin(x+y)', '3*x**4')}
    initial_eta = {"flat": 0.0}
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    V = VectorElement(family='CG', cell=mesh.ufl_cell(), degree=2, dim=2)
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    nmodes = basis_struct.get("nmodes")
    string_V = "["
    string_Q = "["
    for i in range(nmodes-1):
        string_V = string_V + "V,"
        string_Q = string_Q + "Q,"
    string_V = string_V + "V]"
    string_Q = string_Q + "Q]"
    if nmodes == 1:
        W = VectorFunctionSpace(mesh, "CG", 2, dim=2)
        P = FunctionSpace(mesh, "CG", 1)
    else:
        W = FunctionSpace(mesh, MixedElement(eval(string_V)))
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))

    MakeInitialObjectFunc(initial_u, initial_eta, basis_struct, W, P)

    return 1


def test8_degn_vary_flat():

    distname = "uniform"
    sto_poly_deg = 2
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    initial_u = {"vary": ('2*sp.sin(x+y)', '3*x**4')}
    initial_eta = {"vary": '20.*x*y'}
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}

    mesh = MakeMesh(domain)
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    V = VectorElement(family='CG', cell=mesh.ufl_cell(), degree=2, dim=2)
    Q = FiniteElement(family='CG', cell=mesh.ufl_cell(), degree=1)
    nmodes = basis_struct.get("nmodes")
    string_V = "["
    string_Q = "["
    for i in range(nmodes-1):
        string_V = string_V + "V,"
        string_Q = string_Q + "Q,"
    string_V = string_V + "V]"
    string_Q = string_Q + "Q]"
    if nmodes == 1:
        W = VectorFunctionSpace(mesh, "CG", 2, dim=2)
        P = FunctionSpace(mesh, "CG", 1)
    else:
        W = FunctionSpace(mesh, MixedElement(eval(string_V)))
        P = FunctionSpace(mesh, MixedElement(eval(string_Q)))

    MakeInitialObjectFunc(initial_u, initial_eta, basis_struct, W, P)

    return 1


def main():

    print "Don't consider using subclassExpression in IC."

    red = '\033[91m'
    green = '\033[92m'
    blue = '\033[94m'
    bold = '\033[1m'

    flag = test1_deg0_flat_flat()
    if flag:
        print blue+bold+"Test1: "+bold+green+"PASS"
    else:
        print blue+bold+"Test1: "+bold+red+"FAIL"

    flag = test2_deg0_flat_vary()
    if flag:
        print blue+bold+"Test2: "+bold+green+"PASS"
    else:
        print blue+bold+"Test2: "+bold+red+"FAIL"

    flag = test3_deg0_vary_flat()
    if flag:
        print blue+bold+"Test3: "+bold+green+"PASS"
    else:
        print blue+bold+"Test3: "+bold+red+"FAIL"

    flag = test4_deg0_vary_flat()
    if flag:
        print blue+bold+"Test4: "+bold+green+"PASS"
    else:
        print blue+bold+"Test4: "+bold+red+"FAIL"

    flag = test5_degn_flat_flat()
    if flag:
        print blue+bold+"Test5: "+bold+green+"PASS"
    else:
        print blue+bold+"Test5: "+bold+red+"FAIL"

    flag = test6_degn_flat_vary()
    if flag:
        print blue+bold+"Test6: "+bold+green+"PASS"
    else:
        print blue+bold+"Test6: "+bold+red+"FAIL"

    flag = test7_degn_vary_flat()
    if flag:
        print blue+bold+"Test7: "+bold+green+"PASS"
    else:
        print blue+bold+"Test7: "+bold+red+"FAIL"

    flag = test8_degn_vary_flat()
    if flag:
        print blue+bold+"Test8: "+bold+green+"PASS"
    else:
        print blue+bold+"Test8: "+bold+red+"FAIL"

    return 0


if __name__ == '__main__':
    main()
