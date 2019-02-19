from fenics import *
import numpy as np
import chaospy as cp
import sympy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from joblib import Parallel, delayed
import matplotlib.cm as cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import pdb


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
        string = "q"+str(i) + "=sp.Symbol('q"+str(i)+"')"
        exec(string) in globals()

    # create modes one by one.
    length = len(args)
    string = "[sp.integrate( eval(args[k]) * (expr) * JointPDF, "
    for j in range(dim):
        string = string + "(q" + str(j) + ",coefficient[2*"+str(j)+"], coefficient[2*"+str(j)+"+1] )"
        if j < dim-1:
            string = string + ", "
        else:
            string = string + ")"
    string = string + " for expr in eval(str(orthpol))]"
    modes_struct = dict()
    for k in range(length):
        modes = eval(string)
        modes_struct[k] = modes

    return modes_struct

def MakeStoIJK(pc_basis_struct):

    # check the struct contains the expected keys.
    nmodes = pc_basis_struct.get("nmodes", False)
    orthpol = pc_basis_struct.get("basis", False)
    JointCDF = pc_basis_struct.get("jointcdf", False)

    if all([nmodes,orthpol,JointCDF]) is False:
        log(ERROR, "ERROR: obtain NULL value from pc basis struct! Exit program...")
        quit()

    # calculate stoijk based on chaospy.
    stoIJK = np.ndarray(shape=(nmodes,nmodes,nmodes), dtype=float)
    for i in range(nmodes):
        for j in range(nmodes):
            for k in range(nmodes):
                stoIJK[i,j,k] = cp.E(orthpol[i]*orthpol[j]*orthpol[k], JointCDF)
    return stoIJK

# stochastic basis
distname = "uniform"  # "uniform"  or "gaussian"   ==> only enable "uniform" mode at present.
sto_poly_deg = 2 # polynomial chaos order is 2.
sto_poly_dim = 2 # use "q0","q1","q2", ....
coefficient = [0.9, 1.1, 0.8, 1.2] # lower1/upper1--lower2/upper2--...

# stochastic coefficient
sto_viscosity = "2 * x * q1**4 + 3 * y * q1**2 * q0**2"
sto_force = "x * (1 - x) * (4 * x + 12 * y - 3) * q0**4 + y * (1 - y) * (8 * x + 6 * y - 2) * q0 * q1**3"

# boundary condition
sto_boundary = "1.0"


def functioncc(distname, sto_poly_deg, sto_poly_dim, coefficient, sto_viscosity, sto_force, sto_boundary):
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    allmodes = MakeStoModes(basis_struct, sto_viscosity, sto_force)
    v_modes = allmodes.get(0)
    f_modes = allmodes.get(1)
    nmodes = basis_struct.get("nmodes")
    stoIJK = MakeStoIJK(basis_struct)
    boundary_dict = MakeStoModes(basis_struct, sto_boundary)
    bcs_modes = boundary_dict.get(0)
    
    
    mesh = UnitSquareMesh(32, 32)
    V = FiniteElement("CG", mesh.ufl_cell(), 1)
    string_V = "["
    for i in range(nmodes-1):
        string_V = string_V + "V,"
    string_V = string_V + "V]"
    W = FunctionSpace( mesh, MixedElement(eval(string_V)) )
    
    u = TrialFunction(W)
    v = TestFunction(W)
    
    #F_u = inner(grad(u[0]), grad(v[0])) * dx + \
    #    1e-100*inner(grad(u[1]), grad(v[0])) * dx + \
    #    1e-100*inner(grad(u[0]), grad(v[1])) * dx + \
    #    inner(grad(u[1]), grad(v[1])) *dx - \
    #    Constant(1.0)*v[0]*dx - Constant(1.0)*v[1]*dx
    
    F_u = 0.0
    for k in range(nmodes):
        for j in range(nmodes):
            for i in range(nmodes):
                if abs(stoIJK[i][j][k]) >= DOLFIN_EPS:
                    F_u = F_u + stoIJK[i][j][k] * Expression(sp.printing.ccode(v_modes[i]),degree=1) * inner( grad(u[j]), grad(v[k]) ) * dx
        F_u = F_u - Expression(sp.printing.ccode(f_modes[k]), degree=3 ) *v[k]*dx
    
    
    a = assemble(lhs(F_u))
    L = assemble(rhs(F_u))
    #element=W.sub(i).ufl_element()
    def boundaries(x,on_boundary):
        return on_boundary
    #bcs_0 = DirichletBC(W.sub(0), Constant(1.0), boundaries)
    #bcs_1 = DirichletBC(W.sub(1), Constant(2.0), boundaries)
    
    #bcs_expression = [ Expression(sp.printing.ccode(bcs_modes[i]), degree=1)  for i in range(nmodes)]
    bcs = DirichletBC(W, bcs_modes, boundaries)
    bcs.apply(a)
    bcs.apply(L)
    
    
#    u_file = XDMFFile("u.xdmf")
    u0 = Function(W)
    solve(a, u0.vector(), L)
#    u0_list = u0.split(deepcopy=True)
#    CC = FunctionSpace(mesh, "CG", 1)
#    temp = Function(CC, name="u_modes")
#    for i in range(nmodes):
#        temp.assign(u0_list[i])
#        u_file.write(temp, float(Constant(i+1)))

    return u0, basis_struct

nsample = 3
xyticknumber = 6
xytickstep = 4
scatterplotsize = 5

testnodeX = [0.2]
testnodeY = [0.2]
testnodes = [[a,b] for a in testnodeX for b in testnodeY]
samplex = np.linspace(coefficient[0], coefficient[1],  nsample)
sampley = np.linspace(coefficient[2], coefficient[3],  nsample)
sampleX, sampleY = np.meshgrid(samplex, sampley)
u0, basis_struct = functioncc(distname, sto_poly_deg, sto_poly_dim, coefficient, sto_viscosity, sto_force, sto_boundary)
orthpol = basis_struct.get("basis")
nmodes = basis_struct.get("nmodes")

binrandom = np.zeros([nsample,nsample])
for j in range(nsample):
    for k in range(nsample):
        orth_list = [orthpol[mode](samplex[j], sampley[k]) for mode in range(nmodes)]
        binrandom[j,k] = np.dot(orth_list,u0(testnodes[0][0],testnodes[0][1]) )


binrandom_true = np.zeros([nsample,nsample])
for j in range(nsample):
    for k in range(nsample):
        # recal sto_viscosity sto_force
        q0 = samplex[j]
        q1 = sampley[k]
        x=sp.Symbol('x')
        y=sp.Symbol('y')

        sample_viscosity = str(eval(sto_viscosity))
        sample_force = str(eval(sto_force))
        sample_boundary = '1.0'
        u1, basis_struct2 = functioncc(distname, 0, sto_poly_dim, coefficient, sample_viscosity, sample_force, sample_boundary)
        binrandom_true[j,k] = u1(testnodes[0][0], testnodes[0][1])
        print "finish: " + str(j)+str(k)


fig = plt.figure()
ax = fig.gca(projection='3d')
surf = ax.plot_surface(sampleX, sampleY, binrandom, cmap=cm.coolwarm, linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=10)
ax.scatter(sampleX, sampleY, binrandom_true, s=scatterplotsize, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.set_xlim(coefficient[0], coefficient[1])
ax.set_ylim(coefficient[2], coefficient[3])
###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
# ax.set_zlim(1.005,1.013)
ax.set_xticks(np.linspace(coefficient[0], coefficient[1], xyticknumber))
ax.set_yticks(np.linspace(coefficient[2], coefficient[3], xyticknumber))
###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
# ax.set_zticks(np.linspace(1.005,1.013,5))
ax.set_xlabel(r'$\xi_1$')
ax.set_ylabel(r'$\xi_2$')
ax.set_zlabel("u")
plt.title("(" + str(testnodes[0][0]) + ", " + str(testnodes[0][1]) + "); " + "stochastic order=" + str(sto_poly_deg))
# plt.show()
plt.savefig("PCMC_3Dsurface_point_" + str(0) + "_stoorder_" + str(sto_poly_deg) + ".png")
plt.close()
































