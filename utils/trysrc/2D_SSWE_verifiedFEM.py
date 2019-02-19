""" This is a scratch function to help arrange terms in lhs and rhs for this three equations."""


from fenics import *
from mshr import *
import chaospy as cp
import sympy as sp
import numpy as np
import sys
import timeit
import os
from petsc4py import PETSc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from joblib import Parallel, delayed
import matplotlib.cm as cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import pdb

set_log_level(INFO)
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}
sys.setrecursionlimit(10000)
#########################################################
# Function about calculate modes of each expression.
#########################################################
def MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient):

    """
        Function for creating stochastic basis function.
        return: dictionary about basis struct
    """

    # check if the user input is correct.
    if distname == "uniform":
        if len(coefficient)/2 != sto_poly_dim:
            log(ERROR,"ERROR: size of coefficient is not matching the dimension of polynomial chaos basis!")
            quit()
        if any( [ coefficient[2*k] >coefficient[2*k+1] for k in range(len(coefficient)/2) ] ):
            log(ERROR,"ERROR: coefficient value should be in order as lower-high-lower-high, invalid number entering!")
            quit()
    elif distname == "gaussian":
        if len(coefficient)/2 != sto_poly_dim:
            log(ERROR,"ERROR: size of coefficient is not matching the dimension of polynomial chaos basis!")
            quit()
    else:
        log(ERROR, "ERROR: Only implement uniform and gaussian distribution! dist_name not found!")
        quit()

    # start making stochastic basis function.
    if distname == "uniform":
        kercis = [cp.Uniform(coefficient[2*k],coefficient[2*k+1]) for k in range(len(coefficient)/2) ]
        CDFexprString = "cp.J("
        for i in range(len(kercis)):
            CDFexprString = CDFexprString +"kercis["+str(i)+"]"
            if i < len(kercis)-1:
                CDFexprString = CDFexprString + ","
        CDFexprString = CDFexprString + ")"
        JointCDF = eval(CDFexprString)
        JointPDF = 1.0
        for i in range(len(kercis)):
            JointPDF = JointPDF / ( coefficient[2*i+1] - coefficient[2*i] )
        orthpol = cp.orth_ttr(sto_poly_deg, JointCDF, normed=True)
        nmodes= len(orthpol)
    elif distname == "gaussian":
        log(ERROR, "ERROR: Not implement Gaussian distribution yet! Exit program...")
        quit()
    else:
        log(ERROR, "ERROR: Only implement uniform and gaussian distribution! dist_name not found!")
        quit()
    pc_basis_struct = {"name": distname, "dim": sto_poly_dim, "order": sto_poly_deg, "nmodes": nmodes, "coefficient": coefficient, "basis": orthpol, "jointcdf": JointCDF, "jointpdf": JointPDF}

    return pc_basis_struct

def MakeStoModes(pc_basis_struct,*args):
    """
        This function makes the stochastic modes.
    
    :param pc_basis_struct:
    :param args: a set of strings with x, y, q0, q1, q2, ...
    :return: a set of lists with each entry as a string with portable cpp code.
    
    """

    # check the struct contains the expected keys.
    dim = pc_basis_struct.get("dim", False)
    orthpol = pc_basis_struct.get("basis", False)
    JointPDF = pc_basis_struct.get("jointpdf", False)
    coefficient = pc_basis_struct.get("coefficient", False)

    if all([dim,orthpol,JointPDF,coefficient]) is False:
        log(ERROR, "ERROR: obtain NULL value from pc basis struct! Exit program...")
        quit()

    # Initialize chaospy and sympy variables.
    cp.poly.base.POWER = "**"
    cp.poly.base.SEP = "*"
    x=sp.Symbol('x[0]')
    y=sp.Symbol('x[1]')
    t=sp.Symbol('t')
    for i in range(dim):
        str_dim = "q"+str(i) + "=sp.Symbol('q"+str(i)+"')"
        exec(str_dim) in globals()

    # create modes one by one.
    string = "[sp.printing.ccode(  sp.integrate( eval(args[arg]) * (expr) * JointPDF, "
    for j in range(dim):
        string = string + "(q" + str(j) + ",coefficient[2*"+str(j)+"], coefficient[2*"+str(j)+"+1] )"
        if j < dim-1:
            string = string + ", "
        else:
            string = string + ")  )"
    string = string + " for expr in eval(str(orthpol))]"

    # return a single list or nested list
    if len(args) == 1:
        arg = 0
        return eval(string)
    else:
        return [eval(string) for arg in range(len(args))]

def MakeStoIJK(pc_basis_struct):

    # check the struct contains the expected keys.
    nmodes = pc_basis_struct.get("nmodes", False)
    orthpol = pc_basis_struct.get("basis", False)
    JointCDF = pc_basis_struct.get("jointcdf", False)

    if all([nmodes,orthpol,JointCDF]) is False:
        log(ERROR, "ERROR: obtain NULL value from pc basis struct! Exit program...")
        quit()

    # calculate stoijk based on chaospy.
    stoIJK_list = np.ndarray(shape=(nmodes,nmodes,nmodes), dtype=float)
    for i in range(nmodes):
        for j in range(nmodes):
            for k in range(nmodes):
                stoIJK_list[i, j, k] = cp.E(orthpol[i]*orthpol[j]*orthpol[k], JointCDF)
    return stoIJK_list

def MakeMesh(domain_dict):

    # check if the key and the value are entered right.
    key = domain_dict.keys()
    if len(key) != 1:
        log(ERROR, "ERROR: enter more than one domain! not allowed! Exit program...")
        quit()

    # check if the value is as expected.
    value = domain_dict.get(key[0])
    if key[0] == "rectangle":
        if len(value) != 6:
            log(ERROR, "ERROR: expecting 6 domain parameters! Exit program...")
            quit()
    elif key[0] == "importfile":
        if not isinstance(value, str):
            log(ERROR, "ERROR: expecting strings of the mesh full name! Exit program...")
            quit()
    else:
        log(ERROR, "ERROR: not an implemented type of domain! Exit program...")
        quit()

    # Make the mesh.
    if key[0] == "rectangle":
        mesh_object = RectangleMesh(Point(value[0], value[1]), Point(value[2], value[3]), value[4], value[5])
    elif key[0] == "importfile":
        mesh_object = Mesh(value)
    else:
        quit()

    return mesh_object

def MakeBathy(bathy_dict, pc_basis_struct, Bathy_functionspace, current_t):

    # check if the key are entered right.
    key = bathy_dict.keys()
    if len(key) != 1:
        log(ERROR, "ERROR: enter more than one domain! not allowed! Exit program...")
        quit()
    if key[0] != "flat" and key[0] != "vary":
        log(ERROR, "ERROR: don't have this kind of bathymetry implemented! Exit program...")
        quit()

    # check if the value are entered right.
    value = bathy_dict.get(key[0])
    if key[0] == "flat":
        if isinstance(value,float) == False :
            log(ERROR, "ERROR: Wrong bathymetry value entered! Could only take float! Exit program...")
            quit()
    elif key[0] == "vary":
        if ( isinstance(value, str) == False ) and ( value != None ):
            log(ERROR, "ERROR: Wrong bathymetry value entered! Could only take string or 'None'! Exit program...")
            quit()
    else:
        pass

    # make the bathymetry
    if key[0] == "flat":
        bathy_list = MakeStoModes(pc_basis_struct, str(value))
        bathy_function = Constant(bathy_list)
    elif key[0] == "vary":
        if isinstance(value, str):
            x = sp.Symbol('x[0]')
            y = sp.Symbol('x[1]')
            bathy_list = MakeStoModes(pc_basis_struct, value)
            if any(["t" in strlist for strlist in bathy_list]):
                bathy_Expression = Expression(bathy_list, element=Bathy_functionspace.ufl_element(), t=current_t)
            else:
                bathy_Expression = Expression(bathy_list, element=Bathy_functionspace.ufl_element())
        elif value == None:
            log(ERROR,"Not fully implemented yet!! Need to modify the class in order to enable multiple dimension!")
            quit()
            # Need to rewrite the expression in order to make multidimensional values.
            # bathy_Expression = bottomExpression(element=Bathy_functionspace.ufl_element())
        else:
            pass
        bathy_function = interpolate(bathy_Expression, Bathy_functionspace)
    else:
        pass

    return bathy_function

def MakeInitialObjectFunc(initu_dict, initeta_dict, pc_basis_struct, u_functionspace, eta_functionspace):

    # check if the keys is entered right.
    keyu = initu_dict.keys()
    keyeta = initeta_dict.keys()
    if len(keyu) != 1:
        log(ERROR, "ERROR: enter more than one initial condition for u! not allowed! Exit program...")
        quit()
    if len(keyeta) != 1:
        log(ERROR, "ERROR: enter more than one initial condition for eta! not allowed! Exit program...")
        quit()
    if keyu[0] != "flat" and keyu[0] != "vary":
        log(ERROR, "ERROR: don't have this kind of initial condition for u implemented! Exit program...")
        quit()
    if keyeta[0] != "flat" and keyeta[0] != "vary":
        log(ERROR, "ERROR: don't have this kind of initial condition for eta implemented! Exit program...")
        quit()

    # check if the values is entered right
    valueu = initu_dict.get(keyu[0])
    if len(valueu) != 2:
        log(ERROR, "ERROR: enter the wrong length of initial condition for u! Tuple size should be 2! Exit program...")
        quit()
    if isinstance(valueu, tuple) == False:
            log(ERROR, "ERROR: enter the wrong type initial conditions for u! Should be a tuple pair! Exit program...")
            quit()
    if keyu[0] == "flat":
        if any([isinstance(valueu[i], float) == False for i in range(len(valueu))]):
            log(ERROR, "ERROR: enter the wrong type of initial condition for u_vel and v_vel! Should be float for both! Exit program...")
            quit()
    elif keyu[0] == "vary":
        if any([isinstance(valueu[i],str) == False for i in range(len(valueu))]):
            log(ERROR, "ERROR: enter the wrong type initial condition for u! Should be string for both! Exit program...")
            quit()
    else:
        pass
    valueta = initeta_dict.get(keyeta[0])
    if keyeta[0] == "flat":
        if isinstance(valueta, float) == False:
            log(ERROR, "ERROR: enter the wrong type initial condition for eta! Should be a float! Exit program...")
            quit()
    elif keyeta[0] == "vary":
        if isinstance(valueta, str) == False:
            log(ERROR, "ERROR: enter the wrong type initial condition for eta! Should be a string! Exit program...")
            quit()
    else:
        pass

    # make the stochastic initial condition for u and eta.
    nmodes = pc_basis_struct.get("nmodes")
    if keyu[0] == "flat":
        tempu, tempv = MakeStoModes(pc_basis_struct, str(valueu[0]), str(valueu[1]))
    elif keyu[0] == "vary":
        tempu, tempv = MakeStoModes(pc_basis_struct, valueu[0], valueu[1])
    else:
        pass
    if len(tempu) != nmodes or len(tempv) != nmodes:
        log(ERROR, "ERROR: tempu and tempv has wrong length! Should be equal to nmodes! not allowed! Exit program...")
        quit()
    if keyeta[0] == "flat":
        tempeta = MakeStoModes(pc_basis_struct, str(valueta))
    elif keyeta[0] == "vary":
        tempeta = MakeStoModes(pc_basis_struct, valueta)
    else:
        pass
    if len(tempeta) != nmodes:
        log(ERROR, "ERROR: tempeta has wrong length! Should be equal to nmodes! not allowed! Exit program...")
        quit()
    u_ic_list = []
    eta_ic_list = []
    for i in range(nmodes):
        u_ic_list.append(tempu[i])
        u_ic_list.append(tempv[i])
        eta_ic_list.append(tempeta[i])

    u_ic_function = interpolate(Expression(u_ic_list, element=u_functionspace.ufl_element()), u_functionspace)
    if nmodes == 1:
        eta_ic_function = interpolate(Expression(eta_ic_list[0], element=eta_functionspace.ufl_element()), eta_functionspace)
    else:
        eta_ic_function = interpolate(Expression(eta_ic_list, element=eta_functionspace.ufl_element()), eta_functionspace)

    return u_ic_function, eta_ic_function

def MakeBoundaryObjectList(boundary_u_dict, boundary_eta_dict, pc_basis_struct, boundaries, HasFile, u_functionspace, eta_functionspace, current_t):

    # check both u and eta have the correct key==> the key should be subdomain class
    nmodes = pc_basis_struct.get("nmodes")
    nukeys = len(boundary_u_dict.keys())
    netakeys = len(boundary_eta_dict.keys())
    if all([ isinstance(boundary_u_dict.keys()[i], object) for i in range(nukeys) ]) == False:
        log(ERROR, "ERROR: key type of boundary u is wrong! Should be a user-defined class! Exit Program...")
        quit()
    if all([ isinstance(boundary_eta_dict.keys()[i], object) for i in range(netakeys) ]) == False:
        log(ERROR, "ERROR: key type of boundary eta is wrong! Should be a user-defined class! Exit Program...")
        quit()

    # initialize returned bc object list.
    u_bc_list = []
    eta_bc_list = []

    # make boundaries subdomain and mark it sequentially.
    if HasFile:
        pass
    else:
        boundaries.set_all(0)
        allkeys = boundary_u_dict.keys() + boundary_eta_dict.keys()
        for i, key in enumerate(allkeys):
            tempsubdomain = key()
            tempsubdomain.mark(boundaries, i+1)

    # make boundary object list for u.
    for i, key in enumerate(boundary_u_dict.keys()):
        value = boundary_u_dict.get(key)
        if value == "noslip":
            u_list = list(np.zeros(nmodes*2))
            u_bc_list.append(DirichletBC(u_functionspace, u_list, boundaries, i + 1))
        elif value == "freeslipxx":
            if nmodes == 1:
                u_bc_list.append(DirichletBC(u_functionspace.sub(1), Constant(0), boundaries, i + 1))
            else:
                for j in range(nmodes):
                    u_bc_list.append(DirichletBC(u_functionspace.sub(j).sub(1), Constant(0), boundaries, i + 1))
        elif value == "freeslipyy":
            if nmodes == 1:
                u_bc_list.append(DirichletBC(u_functionspace.sub(0), Constant(0), boundaries, i + 1))
            else:
                for j in range(nmodes):
                    u_bc_list.append(DirichletBC(u_functionspace.sub(j).sub(0), Constant(0), boundaries, i + 1))
        elif isinstance(value, tuple) and len(value) == 2 and all([isinstance(value[j], float) for j in range(2)]):
            u_list = list(value) + list(np.zeros((nmodes-1)*2))
            u_bc_list.append(DirichletBC(u_functionspace, u_list, boundaries, i + 1))
        elif isinstance(value, tuple) and len(value) == 2 and all([isinstance(value[j], str) for j in range(2)]):
            tempu_list, tempv_list = MakeStoModes(pc_basis_struct, value[0], value[1])
            u_list = []
            for j in range(nmodes):
                u_list.append(tempu_list[j])
                u_list.append(tempv_list[j])
            if any(["t" in strlist for strlist in u_list]):
                u_list_expression = Expression(u_list, element=u_functionspace.ufl_element(), t=current_t)
                u_bc_list.append(DirichletBC(u_functionspace, u_list_expression, boundaries, i + 1))
            else:
                u_list_expression = Expression(u_list, element=u_functionspace.ufl_element())
                u_bc_list.append(DirichletBC(u_functionspace, u_list_expression, boundaries, i + 1))
        else:
            log(ERROR, "ERROR: enter wrong boundary u type.")
            quit()

    # make boundary object list for eta.
    for i, key in enumerate(boundary_eta_dict.keys()):
        value = boundary_eta_dict.get(key)
        if isinstance(value, float):
            eta_list = [value] + list(np.zeros(nmodes-1))
            if nmodes == 1:
                eta_bc_list.append(DirichletBC(eta_functionspace, eta_list[0], boundaries, i + 1 + nukeys))
            else:
                eta_bc_list.append(DirichletBC(eta_functionspace, eta_list, boundaries, i + 1 + nukeys))
        elif isinstance(value, str):
            eta_list = MakeStoModes(pc_basis_struct, value)
            if any(["t" in strlist for strlist in eta_list]):
                eta_list_expression = Expression(eta_list, element=eta_functionspace.ufl_element(), t=current_t)
                eta_bc_list.append(DirichletBC(eta_functionspace, eta_list_expression, boundaries, i + 1 + nukeys))
            else:
                eta_list_expression = Expression(eta_list, element=eta_functionspace.ufl_element())
                eta_bc_list.append(DirichletBC(eta_functionspace, eta_list_expression, boundaries, i + 1 + nukeys))
        else:
            log(ERROR, "ERROR: enter wrong boundary eta type.")
            quit()

    return u_bc_list, eta_bc_list, u_list, u_list_expression

def test_verify_sfem():

    ##################################################################
    # Input parameter
    ##################################################################

    # stochastic basis
    distname = "uniform"  # "uniform"  or "gaussian"   ==> only enable "uniform" mode at present.
    sto_poly_deg = 2 # polynomial chaos order is 2.
    # the dimension and coefficient should be paired and dim!=0  and coefficient!=Null
    sto_poly_dim = 1 # use "q0","q1","q2", ....
    coefficient = [1,2] # lower1/upper1--lower2/upper2--...

    # horizontal domain setting
    # first way: simple domain built in
    domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
    # second way: complex domain import: domain = {"importfile": "inlet.xml"}

    # stochastic coefficient # if contains sin(), cos(), should use sympy sin ans sympy cos!!!!!!!!!!
    sto_viscosity = "0.0"
    sto_friction = "0.0"
    sto_force_x = "0.0"
    sto_force_y = "0.0"

    # terms control
    include_viscosity = False
    include_advection = False
    linear_divergence = True

    # time parameter setting
    startime = 0.0
    endtime = 25.0
    timestep = 0.5
    theta = 0.5

    # equation parameter setting
    gravity = 9.81

    # bathymetry setting
    bathymetry = {"flat": 20.} # choose between "flat" or "vary"
    #bathymetry = {"vary": "2*x+3*y"}
    #bathymetry = {"vary": None} # when it's "vary", It's either a string "2*x+3*y" or this class should be provided.
    #class bottomExpression(Expression):

    # initial condition
    initial_u = {"flat": (0., 0.)} # choose between "flat" or "vary"
    initial_eta = {"flat": 0.0} # choose between "flat" or "vary"
    #initial_u ={"vary": ('2*x+3*y', '3*y')}
    #initial_eta = {"vary": "20.*sp.sin(x)"}

    # boundary condition
    HasBoundaryfile = False
    # first way: HasBoundaryFile = True
    #boundary_u = {1: ("sp.sin(pi*t/5)*x*(50.0-y)/625.0", "0.0"), 2: "freeslipxx"} # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
    #boundary_eta = {3: 0.0} # "float/str"
    # second way: HasBoundaryFile = False
    class Left(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], domain.get("rectangle")[0])
    class Right(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and near(x[0], domain.get("rectangle")[2])
    class Updown(SubDomain):
        def inside(self, x, on_boundary):
            return on_boundary and (near(x[1], domain.get("rectangle")[1]) or near(x[1], domain.get("rectangle")[3]))
    boundary_u = {Left: ("sp.sin(pi*t/5)*y*(50.0-y)/625.0", "0.0"), Updown: "freeslipxx"}# include sin cos ==> should add sympy symbol!!! # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
    boundary_eta = {Right: 0.0} # "float/str"



    ##################################################################
    # start preparing
    ##################################################################

    # Get stochastic basis
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    stoIJK = MakeStoIJK(basis_struct)

    # Get horizontal domain setting
    mesh = MakeMesh(domain)
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
    B = FunctionSpace(mesh, Q)

    # Get stochastic coefficient
    nu_expression_list, friction_expression_list, f_u_expression_list, f_v_expression_list = MakeStoModes(basis_struct, sto_viscosity, sto_friction, sto_force_x, sto_force_y)

    # Get time parameter setting
    theta = Constant(theta)
    dt = Constant(timestep)
    finish_time = Constant(endtime)
    t = Constant(startime)

    # Get equation parameter settingTrue
    g = Constant(gravity)

    # Get the bathymetry
    project_bottom_function = MakeBathy(bathymetry, basis_struct, P, t)

    # Get initial condition in stochastic setting
    u_ic_function, eta_ic_function = MakeInitialObjectFunc(initial_u, initial_eta, basis_struct, W, P)

    # Get boundary condition in stochastic setting
    if HasBoundaryfile:
        boundaries = MeshFunction("size_t", mesh, "inlet_facet_region.xml")
    else:
        boundaries = FacetFunction("size_t", mesh)
    u_bc_object_list, eta_bc_object_list, u_list, u_list_expression = MakeBoundaryObjectList(boundary_u, boundary_eta, basis_struct, boundaries, HasBoundaryfile, W, P, t)



    ##################################################################
    # starting initializing
    ##################################################################

    # set up solution function
    v = TestFunction(W)
    u = TrialFunction(W)
    q = TestFunction(P)
    eta = TrialFunction(P)

    u00 = Function(W)
    u0 = Function(W, name="u0")
    ut = Function(W) # Tentative velocity
    u1 = Function(W, name="u")
    eta0 = Function(P, name="eta0")
    eta1 = Function(P, name="eta")

    # define total water column
    if linear_divergence:
        H = project_bottom_function
    else:
        H = project_bottom_function + eta0

    # Assign initial condition
    u0.assign(u_ic_function)
    u00.assign(u_ic_function)
    u1.assign(u_ic_function)
    eta0.assign(eta_ic_function)
    eta1.assign(eta_ic_function)

    # For large eddy model ==> error: grad(u) + grad(u).T shape is different, can't match at all. ==> consider: build a new subspace.
    include_les = True

    # Tentative Velocity step
    u_mean = theta * u + (1. - theta) * u0
    u_bash = 3./2 * u0 - 1./2 * u00
    u_diff = u - u0
    F_u_tent = 0.0
    for k in range(nmodes):
        if nmodes == 1:
            F_u_tent = F_u_tent + ( (1/dt) * inner(v, u_diff)* dx +
                                g * inner(v, grad(eta0)) * dx )
        else:
            F_u_tent = F_u_tent + ( (1/dt) * inner(v[2*k], u_diff[2*k])* dx +
                                (1/dt) * inner(v[2*k+1], u_diff[2*k+1]) * dx +
                                g * inner(v[2*k], grad(eta0[k])[0]) * dx +
                                g * inner(v[2*k+1], grad(eta0[k])[1]) * dx )
        for j in range(nmodes):
            for i in range(nmodes):
                if abs(stoIJK[i][j][k]) >= DOLFIN_EPS:
                    F_u_tent = F_u_tent + ( stoIJK[i][j][k] * Expression(friction_expression_list[i], degree=1) * inner(u_mean[2*j], v[2*k]) * dx +
                                            stoIJK[i][j][k] * Expression(friction_expression_list[i], degree=1) * inner(u_mean[2*j+1], v[2*k+1]) * dx )
                    if include_viscosity:
                        F_u_tent = F_u_tent + ( stoIJK[i][j][k] * Expression(nu_expression_list[i], degree=1) * inner(grad(v[2*k]), grad(u_mean[2*j])) * dx +
                                                stoIJK[i][j][k] * Expression(nu_expression_list[i], degree=1) * inner(grad(v[2*k+1]), grad(u_mean[2*j+1])) * dx )
                    if include_advection:
                        F_u_tent = F_u_tent + ( stoIJK[i][j][k] * inner(v[2*k], u_bash[2*i] * grad(u_mean[2*j])[0] + u_bash[2*i+1] * grad(u_mean[2*j])[1]) * dx +
                                                stoIJK[i][j][k] * inner(v[2*k+1], u_bash[2*i] * grad(u_mean[2*j+1])[0] + u_bash[2*i+1] * grad(u_mean[2*j+1])[1]) * dx )
        # if "t" in f_u_expression_list[k]:
        #     tempformu = Expression(f_u_expression_list[k], degree=1, t=Constant(t+theta*dt))
        #     F_u_tent = F_u_tent - v[2*k] * tempformu * dx
        # else:
        #     F_u_tent = F_u_tent - v[2*k] * Expression(f_u_expression_list[k], degree=1) * dx
        # if "t" in f_v_expression_list[k]:
        #     tempformv = Expression(f_v_expression_list[k], degree=1, t=Constant(t+theta*dt))
        #     F_u_tent = F_u_tent - v[2*k+1] * tempformv * dx
        # else:
        #     F_u_tent = F_u_tent - v[2*k+1] * Expression(f_v_expression_list[k], degree=1) * dx
    a_u_tent = lhs(F_u_tent)
    L_u_tent = rhs(F_u_tent)

    # Pressure correction step
    eta_diff = eta - eta0
    ut_mean = theta * ut + (1. - theta) * u0
    F_p_corr = 0.0
    for k in range(nmodes):
        if nmodes == 1:
            F_p_corr = F_p_corr + inner(eta_diff, q) * dx
        else:
            F_p_corr = F_p_corr + inner(eta_diff[k], q[k]) * dx
        for j in range(nmodes):
            for i in range(nmodes):
                if abs(stoIJK[i][j][k]) >= DOLFIN_EPS:
                    if nmodes == 1:
                        F_p_corr = F_p_corr + ( stoIJK[i][j][k] * theta**2 * g * dt**2 * H[0] * inner(grad(q), grad(eta_diff)) * dx \
                                        + stoIJK[i][j][k] * dt * q * grad(H[0]*ut_mean[0])[0] * dx \
                                        + stoIJK[i][j][k] * dt * q * grad(H[0]*ut_mean[1])[1] * dx )
                    else:
                        F_p_corr = F_p_corr + stoIJK[i][j][k] * theta**2 * g * dt**2 * H[i] * inner(grad(q[k]), grad(eta_diff[j])) * dx \
                                        + stoIJK[i][j][k] * dt * q[k] * grad(H[i]*ut_mean[2*j])[0] * dx \
                                        + stoIJK[i][j][k] * dt * q[k] * grad(H[i]*ut_mean[2*j+1])[1] * dx
    a_p_corr = lhs(F_p_corr)
    L_p_corr = rhs(F_p_corr)

    # Velocity correction step
    eta_diff = eta1 - eta0
    F_u_corr = 0.0
    if nmodes == 1:
        F_u_corr = F_u_corr + inner(v, u)*dx - inner(v, ut)*dx + dt*g*theta*inner(v, grad(eta_diff))*dx
    else:
        for k in range(nmodes):
            F_u_corr = F_u_corr + ( u[2*k]*v[2*k] - ut[2*k]*v[2*k] + dt*g*theta*v[2*k]*grad(eta_diff[k])[0] ) * dx \
                            + ( u[2*k+1]*v[2*k+1] - ut[2*k+1]*v[2*k+1] + dt*g*theta*v[2*k+1]*grad(eta_diff[k])[1] ) * dx
    a_u_corr = lhs(F_u_corr)
    L_u_corr = rhs(F_u_corr)

    # Assemble matrices
    A_u_corr = assemble(a_u_corr)
    for bc in u_bc_object_list: bc.apply(A_u_corr)
    a_u_corr_solver = LUSolver(A_u_corr)
    a_u_corr_solver.parameters["reuse_factorization"] = True

    if linear_divergence:
        A_p_corr = assemble(a_p_corr)
        for bc in eta_bc_object_list: bc.apply(A_p_corr)
        a_p_corr_solver = LUSolver(A_p_corr)
        a_p_corr_solver.parameters["reuse_factorization"] = True


    #######################################################################
    # Solving the equation starts ...
    #######################################################################

    # Time stepping
    log(INFO, "Start of time loop")
    timestep = 0
    usinglefs = FunctionSpace(mesh, V)
    etasinglefs = FunctionSpace(mesh, Q)
    tempu1 = Function(usinglefs, name="u_modes")
    tempeta1 = Function(etasinglefs, name="eta_modes")
    #u_file = []
    #eta_file = []

    #u1splits = [0]
    #eta1splits = [0]
    if nmodes == 1:
        #u1splits[0] = u1
        #eta1splits[0] = eta1
        u_file = File("u_0.pvd")
        eta_file = File("eta_0.pvd")
        u_file << (u1, float(t))
        eta_file << (eta1, float(t))
    else:
        u_file = []
        eta_file = []
        u1splits = u1.split(deepcopy=True)
        eta1splits = eta1.split(deepcopy=True)
        for mode in range(nmodes):
            u_file.append(XDMFFile("u_"+str(mode)+".xdmf"))
            eta_file.append(XDMFFile("eta_"+str(mode)+".xdmf"))
            tempu1.assign(u1splits[mode])
            tempeta1.assign(eta1splits[mode])
            u_file[mode].write(tempu1, float(t))
            eta_file[mode].write(tempeta1, float(t))

    while float(t - finish_time) <= - 1e3*DOLFIN_EPS:

        # update current time t.==> successfully updated
        timestep += 1
        t = Constant(t+dt)
        print timestep
        print "time is: ", float(t)

        # update force term with time t. ==> fail to update time.
        if any("t" in strs for strs in f_u_expression_list):
            tempformu.t = Constant(t-(1-theta)*dt)
            print "hehe"
        if any("t" in strs for strs in f_v_expression_list):
            tempformv.t = Constant(t-(1-theta)*dt)
            print "haha"

        # update boundary term with time t. ==> modify time in expression
        u_list_expression.t = t
        #for bcu in u_bc_object_list:
        #    bcu.t = t
        #    print "hehe", float(bcu.t)
        #for bceta in eta_bc_object_list:
        #    bceta.t = t
        #    print "haha", float(bceta.t)

        # update eddy viscosity with time t.


        # Compute tentative velocity step.
        log(INFO, "Solve for tentative velocity.")
        A_u_tent = assemble(a_u_tent)
        b = assemble(L_u_tent)
        for bc in u_bc_object_list: bc.apply(A_u_tent, b)
        solve(A_u_tent, ut.vector(), b)

        # Compute pressure correction step.
        log(INFO, "Solve for pressure correction.")
        b = assemble(L_p_corr)
        for bc in eta_bc_object_list: bc.apply(b)
        if linear_divergence:
            a_p_corr_solver.solve(eta1.vector(), b)
        else:
            A_p_corr = assemble(a_p_corr)
            for bc in eta_bc_object_list: bc.apply(A_p_corr)
            solve(A_p_corr, eta1.vector(), b)

        # Compute Velocity correction step.
        log(INFO, "Solve for velocity update.")
        b = assemble(L_u_corr)
        for bc in u_bc_object_list: bc.apply(b)
        a_u_corr_solver.solve(u1.vector(), b)

        # Rotate functions for next time_step
        u00.assign(u0)
        u0.assign(u1)
        eta0.assign(eta1)

        # Write to file
        if nmodes == 1:
            u_file << (u1, float(t))
            eta_file << (eta1, float(t))
            #u1splits[0] = u1
            #eta1splits[0] = eta1
        else:
            u1splits = u1.split(deepcopy=True)
            eta1splits = eta1.split(deepcopy=True)
            for mode in range(nmodes):
                tempu1.assign(u1splits[mode])
                tempeta1.assign(eta1splits[mode])
                u_file[mode].write(tempu1, float(t))
                eta_file[mode].write(tempeta1, float(t))

        #plot(u1, interactive=True)
        #plot(eta1, interactive=True)

    log(INFO, "End of time loop.")
    print "total time_step is ", timestep

    return 0

if __name__ == "__main__":
    print "This program has been executed, Program Starts ..."
    test_verify_sfem()
else:
    print "This program has been imported ..."
    sys.exit()
