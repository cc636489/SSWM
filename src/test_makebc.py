from fenics import *
from makemesh import MakeMesh
from makestobasis import MakeStoBasis
from makebc import MakeBoundaryObjectList

def test1_deg0_wt_file():

    distname = "uniform"
    t = Constant(10)
    sto_poly_deg = 0
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bcfile = "test1_facet_region.xml"
    if bcfile != None:
        boundary_u = {1: ("sp.sin(pi*t/5)*x*(50.0-y)/625.0*q0", "q1"), 2: "freeslipxx", 3: "freeslipxx"}
        boundary_eta = {4: 0.0}
    else:
        boundary_u = {"Left": ("sp.sin(pi*t/5)*y*(50.0-y)/625.0*q0", "q1"), "Up": "freeslipxx", "Down": "freeslipxx"}# include sin cos ==> should add sympy symbol!!! # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
        boundary_eta = {"Right": 0.0}

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

    u_bc_object_list, eta_bc_object_list, uTimeDependent, etaTimeDependent, u_list_expression, eta_list_expression = MakeBoundaryObjectList(boundary_u, boundary_eta, basis_struct, bcfile, mesh, domain, W, P, t)
    
    print "ut:", uTimeDependent
    print "etat:", etaTimeDependent
    print "ulist:", u_list_expression
    print "etalist:", eta_list_expression
    
    return 1


def test2_deg0_wt_nofile():

    distname = "uniform"
    t = Constant(10)
    sto_poly_deg = 0
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bcfile = None
    if bcfile != None:
        boundary_u = {1: ("sp.sin(pi*t/5)*x*(50.0-y)/625.0*q0", "q1"), 2: "freeslipxx", 3: "freeslipxx"}
        boundary_eta = {4: 0.0}
    else:
        boundary_u = {"Left": ("sp.sin(pi*t/5)*y*(50.0-y)/625.0*q0", "q1"), "Up": "freeslipxx", "Down": "freeslipxx"}# include sin cos ==> should add sympy symbol!!! # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
        boundary_eta = {"Right": 0.0}

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

    u_bc_object_list, eta_bc_object_list, uTimeDependent, etaTimeDependent, u_list_expression, eta_list_expression = MakeBoundaryObjectList(boundary_u, boundary_eta, basis_struct, bcfile, mesh, domain, W, P, t)
    
    print "ut:", uTimeDependent
    print "etat:", etaTimeDependent
    print "ulist:", u_list_expression
    print "etalist:", eta_list_expression

    return 1


def test3_deg0_wot_file():

    distname = "uniform"
    t = Constant(10)
    sto_poly_deg = 0
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bcfile = "test1_facet_region.xml"
    if bcfile != None:
        boundary_u = {1: ("sp.sin(pi/5)*x*(50.0-y)/625.0*q0", "q1"), 2: "freeslipxx", 3: "freeslipxx"}
        boundary_eta = {4: 0.0}
    else:
        boundary_u = {"Left": ("sp.sin(pi/5)*y*(50.0-y)/625.0*q0", "q1"), "Up": "freeslipxx", "Down": "freeslipxx"}# include sin cos ==> should add sympy symbol!!! # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
        boundary_eta = {"Right": 0.0}

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

    u_bc_object_list, eta_bc_object_list, uTimeDependent, etaTimeDependent, u_list_expression, eta_list_expression = MakeBoundaryObjectList(boundary_u, boundary_eta, basis_struct, bcfile, mesh, domain, W, P, t)
    
    print "ut:", uTimeDependent
    print "etat:", etaTimeDependent
    print "ulist:", u_list_expression
    print "etalist:", eta_list_expression

    return 1


def test4_deg0_wot_nofile():

    distname = "uniform"
    t = Constant(10)
    sto_poly_deg = 0
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bcfile = None
    if bcfile != None:
        boundary_u = {1: ("sp.sin(pi/5)*x*(50.0-y)/625.0*q0", "q1"), 2: "freeslipxx", 3: "freeslipxx"}
        boundary_eta = {4: 0.0}
    else:
        boundary_u = {"Left": ("sp.sin(pi/5)*y*(50.0-y)/625.0*q0", "q1"), "Up": "freeslipxx", "Down": "freeslipxx"}# include sin cos ==> should add sympy symbol!!! # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
        boundary_eta = {"Right": 0.0}

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

    u_bc_object_list, eta_bc_object_list, uTimeDependent, etaTimeDependent, u_list_expression, eta_list_expression = MakeBoundaryObjectList(boundary_u, boundary_eta, basis_struct, bcfile, mesh, domain, W, P, t)
    
    print "ut:", uTimeDependent
    print "etat:", etaTimeDependent
    print "ulist:", u_list_expression
    print "etalist:", eta_list_expression

    return 1


def test5_degn_wt_file():

    distname = "uniform"
    t = Constant(10)
    sto_poly_deg = 2
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bcfile = "test1_facet_region.xml"
    if bcfile != None:
        boundary_u = {1: ("sp.sin(pi*t/5)*x*(50.0-y)/625.0*q0", "q1"), 2: "freeslipxx", 3: "freeslipxx"}
        boundary_eta = {4: 0.0}
    else:
        boundary_u = {"Left": ("sp.sin(pi*t/5)*y*(50.0-y)/625.0*q0", "q1"), "Up": "freeslipxx", "Down": "freeslipxx"}# include sin cos ==> should add sympy symbol!!! # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
        boundary_eta = {"Right": 0.0}

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

    u_bc_object_list, eta_bc_object_list, uTimeDependent, etaTimeDependent, u_list_expression, eta_list_expression = MakeBoundaryObjectList(boundary_u, boundary_eta, basis_struct, bcfile, mesh, domain, W, P, t)
    
    print "ut:", uTimeDependent
    print "etat:", etaTimeDependent
    print "ulist:", u_list_expression
    print "etalist:", eta_list_expression
    
    return 1


def test6_degn_wt_nofile():

    distname = "uniform"
    t = Constant(10)
    sto_poly_deg = 2
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bcfile = None
    if bcfile != None:
        boundary_u = {1: ("sp.sin(pi*t/5)*x*(50.0-y)/625.0*q0", "q1"), 2: "freeslipxx", 3: "freeslipxx"}
        boundary_eta = {4: 0.0}
    else:
        boundary_u = {"Left": ("sp.sin(pi*t/5)*y*(50.0-y)/625.0*q0", "q1"), "Up": "freeslipxx", "Down": "freeslipxx"}# include sin cos ==> should add sympy symbol!!! # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
        boundary_eta = {"Right": 0.0}

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

    u_bc_object_list, eta_bc_object_list, uTimeDependent, etaTimeDependent, u_list_expression, eta_list_expression = MakeBoundaryObjectList(boundary_u, boundary_eta, basis_struct, bcfile, mesh, domain, W, P, t)
    
    print "ut:", uTimeDependent
    print "etat:", etaTimeDependent
    print "ulist:", u_list_expression
    print "etalist:", eta_list_expression

    return 1


def test7_degn_wot_file():

    distname = "uniform"
    t = Constant(10)
    sto_poly_deg = 2
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bcfile = "test1_facet_region.xml"
    if bcfile != None:
        boundary_u = {1: ("sp.sin(pi/5)*x*(50.0-y)/625.0*q0", "q1"), 2: "freeslipxx", 3: "freeslipxx"}
        boundary_eta = {4: 0.0}
    else:
        boundary_u = {"Left": ("sp.sin(pi/5)*y*(50.0-y)/625.0*q0", "q1"), "Up": "freeslipxx", "Down": "freeslipxx"}# include sin cos ==> should add sympy symbol!!! # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
        boundary_eta = {"Right": 0.0}

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

    u_bc_object_list, eta_bc_object_list, uTimeDependent, etaTimeDependent, u_list_expression, eta_list_expression = MakeBoundaryObjectList(boundary_u, boundary_eta, basis_struct, bcfile, mesh, domain, W, P, t)
    
    print "ut:", uTimeDependent
    print "etat:", etaTimeDependent
    print "ulist:", u_list_expression
    print "etalist:", eta_list_expression

    return 1


def test8_degn_wot_nofile():

    distname = "uniform"
    t = Constant(10)
    sto_poly_deg = 2
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    bcfile = None
    if bcfile != None:
        boundary_u = {1: ("sp.sin(pi/5)*x*(50.0-y)/625.0*q0", "q1"), 2: "freeslipxx", 3: "freeslipxx"}
        boundary_eta = {4: 0.0}
    else:
        boundary_u = {"Left": ("sp.sin(pi/5)*y*(50.0-y)/625.0*q0", "q1"), "Up": "freeslipxx", "Down": "freeslipxx"}# include sin cos ==> should add sympy symbol!!! # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
        boundary_eta = {"Right": 0.0}

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

    u_bc_object_list, eta_bc_object_list, uTimeDependent, etaTimeDependent, u_list_expression, eta_list_expression = MakeBoundaryObjectList(boundary_u, boundary_eta, basis_struct, bcfile, mesh, domain, W, P, t)
    
    print "ut:", uTimeDependent
    print "etat:", etaTimeDependent
    print "ulist:", u_list_expression
    print "etalist:", eta_list_expression

    return 1


def main():

    print "Don't consider using subclassExpression in BC."

    red = '\033[91m'
    green = '\033[92m'
    blue = '\033[94m'
    bold = '\033[1m'

    flag = test1_deg0_wt_file()
    if flag:
        print blue+bold+"Test1: "+bold+green+"PASS"
    else:
        print blue+bold+"Test1: "+bold+red+"FAIL"

    flag = test2_deg0_wt_nofile()
    if flag:
        print blue+bold+"Test2: "+bold+green+"PASS"
    else:
        print blue+bold+"Test2: "+bold+red+"FAIL"

    flag = test3_deg0_wot_file()
    if flag:
        print blue+bold+"Test3: "+bold+green+"PASS"
    else:
        print blue+bold+"Test3: "+bold+red+"FAIL"

    flag = test4_deg0_wot_nofile()
    if flag:
        print blue+bold+"Test4: "+bold+green+"PASS"
    else:
        print blue+bold+"Test4: "+bold+red+"FAIL"

    flag = test5_degn_wt_file()
    if flag:
        print blue+bold+"Test5: "+bold+green+"PASS"
    else:
        print blue+bold+"Test5: "+bold+red+"FAIL"

    flag = test6_degn_wt_nofile()
    if flag:
        print blue+bold+"Test6: "+bold+green+"PASS"
    else:
        print blue+bold+"Test6: "+bold+red+"FAIL"

    flag = test7_degn_wot_file()
    if flag:
        print blue+bold+"Test7: "+bold+green+"PASS"
    else:
        print blue+bold+"Test7: "+bold+red+"FAIL"

    flag = test8_degn_wot_nofile()
    if flag:
        print blue+bold+"Test8: "+bold+green+"PASS"
    else:
        print blue+bold+"Test8: "+bold+red+"FAIL"

    return 0


if __name__ == '__main__':
    main()
