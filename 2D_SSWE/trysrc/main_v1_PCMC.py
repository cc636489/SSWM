from fenics import *
import sys
import sympy as sp
import numpy as np
from input_test4 import *
from makestobasis import MakeStoBasis
from makestomodes import MakeStoModes
from makestoijk import MakeStoIJK
from makemesh import MakeMesh
from make_bath import make_bath
from makeic import MakeInitialObjectFunc
from makebc import MakeBoundaryObjectList
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import matplotlib.cm as cm


#############################################################
#############################################################
#############################################################
# there is a bug in MC method calculation!! CHECKINg debug.py
#############################################################
#############################################################
#############################################################


set_log_level(INFO)
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True,
               "eliminate_zeros": True,
               "precompute_basis_const": True,
               "precompute_ip_const": True}
sys.setrecursionlimit(10000)

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

# Get stochastic coefficient
nu_expression_list, friction_expression_list, f_u_expression_list, f_v_expression_list = MakeStoModes(basis_struct, sto_viscosity, sto_friction, sto_force_x, sto_force_y)

# Get time parameter setting
theta = Constant(theta)
dt = Constant(timestep)
finish_time = Constant(endtime)
t = Constant(startime)
ntimestep = int( (endtime-startime) / timestep )

# Get equation parameter settingTrue
g = Constant(gravity)

# Get the bathymetry
project_bottom_function = make_bath(bathymetry, basis_struct, P)

# Get initial condition in stochastic setting
u_ic_function, eta_ic_function = MakeInitialObjectFunc(initial_u, initial_eta, basis_struct, W, P)

# Get boundary condition in stochastic setting
u_bc_object_list, eta_bc_object_list, uTimeDependent, etaTimeDependent, u_list_expression, eta_list_expression = MakeBoundaryObjectList(boundary_u, boundary_eta, basis_struct, bcfile, mesh, domain, W, P, t, project_bottom_function)


##################################################################
# starting initializing
##################################################################

# For large eddy model ==> error: grad(u) + grad(u).T shape is different, can't match at all. ==> consider: build a new subspace.
class LES(object):
    
    def __init__(self, V, u, smagorinsky_coefficient):
        self.V = V
        self.u = u
        self.smagorinsky_coefficient = smagorinsky_coefficient
        self.eddy_viscosity = Function(V)
    
    def _strain_rate_tensor(self):
        S = 0.5 * (grad(self.u) + grad(self.u).T)
        return S
    
    def _eddy_viscosity_eqn(self):
        
        dim = len(self.u)
        w = TestFunction(self.V)
        q = TrialFunction(self.V)
        
        cell_vol = CellVolume(self.V.mesh())
        filter_width = cell_vol ** (1.0 / dim)
        
        S = self._strain_rate_tensor()
        second_invariant = 0.0
        for i in range(0, dim):
            for j in range(0, dim):
                second_invariant += 2.0 * (S[i, j] ** 2)
        
        second_invariant = sqrt(second_invariant)
        form = (self.smagorinsky_coefficient * filter_width) ** 2 * second_invariant
        
        lhs = inner(w, q) * dx
        rhs = inner(w, form) * dx
        
        return lhs, rhs
    
    def solve(self):
        # Create a eddy viscosity solver
        les_lhs, les_rhs = self._eddy_viscosity_eqn()
        
        eddy_viscosity_problem = LinearVariationalProblem(les_lhs, les_rhs,
                                                          self.eddy_viscosity, bcs=[])
        _solver = LinearVariationalSolver(eddy_viscosity_problem)
        _solver.parameters["linear_solver"] = "lu"
        _solver.parameters["symmetric"] = True
        _solver.parameters["lu_solver"]["reuse_factorization"] = True
        _solver.solve()
        return self.eddy_viscosity

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

# set up the model parameters in stochastic space
A = FunctionSpace(mesh, V)
B = FunctionSpace(mesh, Q)
f_u_expression_object_list = []
f_v_expression_object_list = []
friction_expression_object_list = []
nu_expression_object_list = []
les = []
for k in range(nmodes):
    if "t" in f_u_expression_list[k]:
        f_u_expression_object_list.append(Expression(f_u_expression_list[k], element=A.sub(0).ufl_element(), t=t))
    else:
        f_u_expression_object_list.append(Expression(f_u_expression_list[k], element=A.sub(0).ufl_element()))
    if "t" in f_v_expression_list[k]:
        f_v_expression_object_list.append(Expression(f_v_expression_list[k], element=A.sub(1).ufl_element(), t=t))
    else:
        f_v_expression_object_list.append(Expression(f_v_expression_list[k], element=A.sub(1).ufl_element()))
    friction_expression_object_list.append(Expression(friction_expression_list[k], element=B.ufl_element()))

    if include_les:
        if nmodes == 1:
            les.append(LES(B, u0, les_parameters['smagorinsky_coefficient']))
        else:
            les.append(LES(B, u0.split()[k], les_parameters['smagorinsky_coefficient']))
        nu_expression_object_list.append(Expression(nu_expression_list[k], element=B.ufl_element()) + les[k].eddy_viscosity)

        #if DEBUG_mode:
            #p = plot(project(les[k].eddy_viscosity, B), title='synthesize viscousity')
            #p.set_cmap('jet')
            #plt.colorbar(p)
            #plt.close()

    else:
        nu_expression_object_list.append(Expression(nu_expression_list[k], element=B.ufl_element()))


# Tentative Velocity step
u_mean = theta * u + (1. - theta) * u0
u_bash = 3./2 * u0 - 1./2 * u00
u_diff = u - u0
F_u_tent = 0.0


for k in range(nmodes):
    if nmodes == 1:
        F_u_tent = F_u_tent + (1/dt) * inner(v, u_diff)*dx +\
            g * inner(v, grad(eta0)) * dx
    else:
        F_u_tent = F_u_tent + (1/dt) * inner(v[2*k], u_diff[2*k]) * dx +\
                            (1/dt) * inner(v[2*k+1], u_diff[2*k+1]) * dx +\
                            g * inner(v[2*k], grad(eta0[k])[0]) * dx +\
                            g * inner(v[2*k+1], grad(eta0[k])[1]) * dx
    for j in range(nmodes):
        for i in range(nmodes):
            if abs(stoIJK[i][j][k]) > 1e-6:
                F_u_tent = F_u_tent + (stoIJK[i][j][k] * friction_expression_object_list[i] * inner(u_mean[2*j], v[2*k]) * dx +
                                        stoIJK[i][j][k] * friction_expression_object_list[i] * inner(u_mean[2*j+1], v[2*k+1]) * dx)
                if include_viscosity:
                    F_u_tent = F_u_tent + (stoIJK[i][j][k] * nu_expression_object_list[i] * inner(grad(v[2*k]), grad(u_mean[2*j])) * dx +
                                            stoIJK[i][j][k] * nu_expression_object_list[i] * inner(grad(v[2*k+1]), grad(u_mean[2*j+1])) * dx)
                if include_advection:
                    F_u_tent = F_u_tent + ( stoIJK[i][j][k] * inner(v[2*k], u_bash[2*i] * grad(u_mean[2*j])[0] + u_bash[2*i+1] * grad(u_mean[2*j])[1]) * dx +
                                            stoIJK[i][j][k] * inner(v[2*k+1], u_bash[2*i] * grad(u_mean[2*j+1])[0] + u_bash[2*i+1] * grad(u_mean[2*j+1])[1]) * dx )
    F_u_tent = F_u_tent - v[2*k] * f_u_expression_object_list[k] * dx
    F_u_tent = F_u_tent - v[2*k+1] * f_v_expression_object_list[k] * dx
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
                    F_p_corr = F_p_corr + (stoIJK[i][j][k] * theta**2 * g * dt**2 * H * inner(grad(q), grad(eta_diff)) * dx
                                    + stoIJK[i][j][k] * dt * q * grad(H*ut_mean[0])[0] * dx
                                    + stoIJK[i][j][k] * dt * q * grad(H*ut_mean[1])[1] * dx)
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
    F_u_corr = inner(v, u)*dx - inner(v, ut)*dx + dt*g*theta*inner(v, grad(eta_diff))*dx
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
timestepcount = 0
usinglefs = FunctionSpace(mesh, V)
etasinglefs = FunctionSpace(mesh, Q)
tempu1 = Function(usinglefs, name="u_modes")
tempeta1 = Function(etasinglefs, name="eta_modes")

if nmodes == 1:
    u_file = File(outputdir + "u_" + outputstr + "only_mode_0.pvd")
    eta_file = File(outputdir + "eta_" + outputstr + "only_mode_0.pvd")
    eddy_file = File(outputdir + "eddy_" + outputstr + "only_mode_0.pvd")
    u_file << (u1, float(t))
    eta_file << (eta1, float(t))
else:
    u_file = []
    eta_file = []
    u1splits = u1.split(deepcopy=True)
    eta1splits = eta1.split(deepcopy=True)
    for mode in range(nmodes):

        if USEpvd:
            u_file.append(File(outputdir + "u_" + outputstr + "{:02d}".format(mode) + ".pvd"))
            eta_file.append(File(outputdir + "eta_" + outputstr + "{:02d}".format(mode) + ".pvd"))
        else:
            u_file.append(XDMFFile(outputdir + "u_"+outputstr+"{:02d}".format(mode)+".xdmf"))
            eta_file.append(XDMFFile(outputdir + "eta_"+outputstr+"{:02d}".format(mode)+".xdmf"))

        tempu1.assign(u1splits[mode])
        tempeta1.assign(eta1splits[mode])

        if USEpvd:
            u_file[mode] << (tempu1, float(t))
            eta_file[mode] << (tempeta1, float(t))
        else:
            u_file[mode].write(tempu1, float(t))
            eta_file[mode].write(tempeta1, float(t))

# Prepare binrandom for PC
testnodes = [[a,b] for a in testnodeX for b in testnodeY]
samplex = np.linspace(coefficient[0], coefficient[1],  nsample)
sampley = np.linspace(coefficient[2], coefficient[3],  nsample)
sampleX, sampleY = np.meshgrid(samplex, sampley)
binrandom_u1 = np.zeros([nsample, nsample, ntimestep+1])
binrandom_v1 = np.zeros([nsample, nsample, ntimestep+1])
binrandom_eta1 = np.zeros([nsample, nsample, ntimestep+1])
orthpol = basis_struct.get("basis")
u1_list = [u1(testnodes[0][0], testnodes[0][1])[2*j] for j in range(nmodes)]
v1_list = [u1(testnodes[0][0], testnodes[0][1])[2*j+1] for j in range(nmodes)]
for j in range(nsample):
    for k in range(nsample):
        orth_list = [orthpol[mode](samplex[j], sampley[k]) for mode in range(nmodes)]
        binrandom_u1[j, k, timestepcount] = np.dot(orth_list, u1_list)
        binrandom_v1[j, k, timestepcount] = np.dot(orth_list, v1_list)
        binrandom_eta1[j, k, timestepcount] = np.dot(orth_list, eta1(testnodes[0][0], testnodes[0][1]))

# Start timestepping.
while float(t - finish_time) < -1e-3:

    # update current time t.==> successfully updated
    timestepcount += 1
    t = Constant(t+dt)
    print timestepcount
    print "time is: ", float(t)

    # update force term with time t. ==> successfully updated. ==> weak force implemented here.
    t_theta = Constant(t-(1-theta)*dt)
    for k in range(nmodes):
        f_u_expression_object_list[k].t = t_theta
        f_v_expression_object_list[k].t = t_theta

    # update boundary term with time t. ==> successfully updated. ==> strong dirichlet BC implemented here.
    if uTimeDependent:
        u_list_expression.t = t
    if etaTimeDependent:
        eta_list_expression.t = t

    # update eddy_viscosity
    if include_les:
        log(INFO, "Compute eddy viscosity.")
        for k in range(nmodes):
            if nmodes == 1:
                les[k].u = u0
            else:
                les[k].u = u0.split()[k]
            les[k].solve()
            #nu_expression_object_list[k] = Expression(nu_expression_list[k], element=B.ufl_element()) + les[k].eddy_viscosity

            if DEBUG_mode:
                cc = project(nu_expression_object_list[k], B)
                for i in range(800, 1010, 10):
                    for j in range(180, 210, 10):
                        print("i = ", i, "j = ", j, cc(i, j))


   # Compute tentative velocity step.
    log(INFO, "Solve for tentative velocity.")
    A_u_tent = assemble(a_u_tent)
    
    if DEBUG_mode:
        print("utent matrix:")
        print(A_u_tent.getrow(0))
    
    b = assemble(L_u_tent)
    
    if DEBUG_mode:
        print("utent vector:")
        print(b.array())
    
    for bc in u_bc_object_list: bc.apply(A_u_tent, b)
    solve(A_u_tent, ut.vector(), b)



    # Compute pressure correction step.
    log(INFO, "Solve for pressure correction.")
    b = assemble(L_p_corr)
    
    if DEBUG_mode:
        print("pcorr vector:")
        print(b.array())
    
    for bc in eta_bc_object_list: bc.apply(b)
    if linear_divergence:
        a_p_corr_solver.solve(eta1.vector(), b)
    else:
        A_p_corr = assemble(a_p_corr)
        
        if DEBUG_mode:
            print("pcorr matrix:")
            print(A_p_corr.getrow(0))
        
        for bc in eta_bc_object_list: bc.apply(A_p_corr)
        solve(A_p_corr, eta1.vector(), b)



    # Compute Velocity correction step.
    log(INFO, "Solve for velocity update.")
    b = assemble(L_u_corr)
    
    if DEBUG_mode:
        print("ucorr vector:")
        print(b.array()[-10:-1])
    
    for bc in u_bc_object_list: bc.apply(b)
    
    if DEBUG_mode:
        print("ucorr matrix:")
        print(A_u_corr.getrow(0))
    
    a_u_corr_solver.solve(u1.vector(), b)



    # Rotate functions for next time_step
    u00.assign(u0)
    u0.assign(u1)
    eta0.assign(eta1)
    
    # if DEBUG_mode:
    #    for i in range(500, 1010, 10):
    #        for j in range(0, 210, 10):
    #            print("i = ", i, "j = ", j, "u1=", u1(i, j))
    #    for i in range(500, 1010, 10):
    #         for j in range(0, 210, 10):
    #             print("i = ", i, "j = ", j, "eta1 = ", eta1(i, j))

    print("===========================================")



    # Write to file
    if nmodes == 1:
        u_file << (u1, float(t))
        eta_file << (eta1, float(t))
    else:
        u1splits = u1.split(deepcopy=True)
        eta1splits = eta1.split(deepcopy=True)
        for mode in range(nmodes):
            
            tempu1.assign(u1splits[mode])
            tempeta1.assign(eta1splits[mode])
            
            if USEpvd:
                u_file[mode] << (tempu1, float(t))
                eta_file[mode] << (tempeta1, float(t))
            else:
                u_file[mode].write(tempu1, float(t))
                eta_file[mode].write(tempeta1, float(t))

    # Calculate bin_random_u1 and bin_random_eta1
    u1_list = [u1(testnodes[0][0], testnodes[0][1])[2*j] for j in range(nmodes)]
    v1_list = [u1(testnodes[0][0], testnodes[0][1])[2*j+1] for j in range(nmodes)]
    for j in range(nsample):
        for k in range(nsample):
            orth_list = [orthpol[mode](samplex[j], sampley[k]) for mode in range(nmodes)]
            binrandom_u1[j, k, timestepcount] = np.dot(orth_list, u1_list)
            binrandom_v1[j, k, timestepcount] = np.dot(orth_list, v1_list)
            binrandom_eta1[j, k, timestepcount] = np.dot(orth_list, eta1(testnodes[0][0], testnodes[0][1]))
            print "PC finish: " + str(j)+str(k)

log(INFO, "End of time loop.")
print "total time_step is ", timestepcount


# Start calculating MC( order = 0 ) here. Will use global variables. but will not change it inside for sure.
def func_mc(sample_vs, sample_fric, sample_f_x, sample_f_y, sample_bathymetry, sample_ic_u, sample_ic_eta, sample_bc_u, sample_bc_eta):

    t_mc = Constant(startime)
    basis_struct_mc = MakeStoBasis(distname, 0, sto_poly_dim, coefficient)

    # Get horizontal domain setting
    W_mc = VectorFunctionSpace(mesh, "CG", 2, dim=2)
    P_mc = FunctionSpace(mesh, "CG", 1)

    # Get stochastic coefficient
    nu_expression_list_mc, friction_expression_list_mc, f_u_expression_list_mc, f_v_expression_list_mc = \
        MakeStoModes(basis_struct_mc, sample_vs, sample_fric, sample_f_x, sample_f_y)
    
    # Get the bathymetry
    project_bottom_function_mc = make_bath(sample_bathymetry, basis_struct_mc, P_mc)
    
    # Get initial condition in stochastic setting
    u_ic_function_mc, eta_ic_function_mc = MakeInitialObjectFunc(sample_ic_u, sample_ic_eta, basis_struct_mc, W_mc, P_mc)
    
    # Get boundary condition in stochastic setting
    u_bc_object_list_mc, eta_bc_object_list_mc, uTimeDependent_mc, etaTimeDependent_mc, u_list_expression_mc, eta_list_expression_mc = MakeBoundaryObjectList(sample_bc_u, sample_bc_eta, basis_struct_mc, bcfile, mesh, domain, W_mc, P_mc, t_mc, project_bottom_function_mc)
    
    ##################################################################
    # starting initializing
    ##################################################################

    # set up solution function
    v_mc = TestFunction(W_mc)
    u_mc = TrialFunction(W_mc)
    q_mc = TestFunction(P_mc)
    eta_mc = TrialFunction(P_mc)
    
    u00_mc = Function(W_mc)
    u0_mc = Function(W_mc, name="u0")
    ut_mc = Function(W_mc) # Tentative velocity
    u1_mc = Function(W_mc, name="u")
    eta0_mc = Function(P_mc, name="eta0")
    eta1_mc = Function(P_mc, name="eta")

    # define total water column
    if linear_divergence:
        H_mc = project_bottom_function_mc
    else:
        H_mc = project_bottom_function_mc + eta0_mc
    
    # Assign initial condition
    u0_mc.assign(u_ic_function_mc)
    u00_mc.assign(u_ic_function_mc)
    u1_mc.assign(u_ic_function_mc)
    eta0_mc.assign(eta_ic_function_mc)
    eta1_mc.assign(eta_ic_function_mc)

    # about fu fv and nu
    if "t" in f_u_expression_list_mc[0]:
        f_u_expression_object_list_mc = Expression(f_u_expression_list_mc[0], element=W_mc.sub(0).ufl_element(), t=t_mc)
    else:
        f_u_expression_object_list_mc = Expression(f_u_expression_list_mc[0], element=W_mc.sub(0).ufl_element())
    if "t" in f_v_expression_list_mc[0]:
        f_v_expression_object_list_mc = Expression(f_v_expression_list_mc[0], element=W_mc.sub(1).ufl_element(), t=t_mc)
    else:
        f_v_expression_object_list_mc = Expression(f_v_expression_list_mc[0], element=W_mc.sub(1).ufl_element())
    friction_expression_object_list_mc = Expression(friction_expression_list_mc[0], element=P_mc.ufl_element())
    nu_expression_mc = Expression(nu_expression_list_mc[0], element=P_mc.ufl_element())

    if include_les:
        les_V_mc = FunctionSpace(mesh, "CG", 1)
        les_mc = LES(les_V_mc, u0_mc, les_parameters['smagorinsky_coefficient'])
        nu_mc = nu_expression_mc + les_mc.eddy_viscosity
    else:
        nu_mc = nu_expression_mc


    # Tentative Velocity step
    u_mean_mc = theta * u_mc + (1. - theta) * u0_mc
    u_bash_mc = 3./2 * u0_mc - 1./2 * u00_mc
    u_diff_mc = u_mc - u0_mc


    F_u_tent_mc = (1/dt) * inner(v_mc, u_diff_mc) * dx + g * inner(v_mc, grad(eta0_mc)) * dx \
        + friction_expression_object_list_mc * inner(u_mean_mc, v_mc) * dx \
        - v_mc[0] * f_u_expression_object_list_mc * dx - v_mc[1] * f_v_expression_object_list_mc * dx

    if include_viscosity:
        F_u_tent_mc += nu_mc * inner(grad(v_mc), grad(u_mean_mc)) * dx

    if include_advection:
        F_u_tent_mc += inner(v_mc, grad(u_mean_mc)*u_bash_mc) * dx

    a_u_tent_mc = lhs(F_u_tent_mc)
    L_u_tent_mc = rhs(F_u_tent_mc)



    # Pressure correction step
    eta_diff_mc = eta_mc - eta0_mc
    ut_mean_mc = theta * ut_mc + (1. - theta) * u0_mc
    F_p_corr_mc = eta_diff_mc*q_mc * dx + theta**2 * g * dt**2 * H_mc * inner(grad(q_mc), grad(eta_diff_mc)) * dx \
        + dt * q_mc * div(H_mc*ut_mean_mc) * dx
    a_p_corr_mc = lhs(F_p_corr_mc)
    L_p_corr_mc = rhs(F_p_corr_mc)



    # Velocity correction step
    eta_diff_mc = eta1_mc - eta0_mc
    F_u_corr_mc = inner(v_mc, u_mc)*dx - inner(v_mc, ut_mc)*dx + dt*g*theta*inner(v_mc, grad(eta_diff_mc))*dx
    a_u_corr_mc = lhs(F_u_corr_mc)
    L_u_corr_mc = rhs(F_u_corr_mc)



    # Assemble matrices
    A_u_corr_mc = assemble(a_u_corr_mc)
    for bc_mc in u_bc_object_list_mc: bc_mc.apply(A_u_corr_mc)
    a_u_corr_solver_mc = LUSolver(A_u_corr_mc)
    a_u_corr_solver_mc.parameters["reuse_factorization"] = True
    
    if linear_divergence:
        A_p_corr_mc = assemble(a_p_corr_mc)
        for bc_mc in eta_bc_object_list_mc: bc_mc.apply(A_p_corr_mc)
        a_p_corr_solver_mc = LUSolver(A_p_corr_mc)
        a_p_corr_solver_mc.parameters["reuse_factorization"] = True

    #######################################################################
    # Solving the equation starts ...
    #######################################################################

    # Time stepping
    timestepcount_mc = 0
    result_u1 = []
    result_v1 = []
    result_eta1 = []
    result_u1.append(u1_mc(testnodes[0][0], testnodes[0][1])[0])
    result_v1.append(u1_mc(testnodes[0][0], testnodes[0][1])[1])
    result_eta1.append(eta1_mc(testnodes[0][0], testnodes[0][1]))

    while float(t_mc - finish_time) < -1e-3:

        # update current time t.==> successfully updated
        timestepcount_mc += 1
        t_mc = Constant(t_mc+dt)

        # update force term with time t. ==> successfully updated. ==> weak force implemented here.
        t_theta_mc = Constant(t_mc-(1-theta)*dt)
        f_u_expression_object_list_mc.t = t_theta_mc
        f_v_expression_object_list_mc.t = t_theta_mc

        # update boundary term with time t. ==> successfully updated. ==> strong dirichlet BC implemented here.
        if uTimeDependent_mc:
            u_list_expression_mc.t = t_mc
        if etaTimeDependent_mc:
            eta_list_expression_mc.t = t_mc

        # update eddy_viscosity
        if include_les:
            les_mc.u = u0_mc
            les_mc.solve()

        # Compute tentative velocity step.
        A_u_tent_mc = assemble(a_u_tent_mc)
        b_mc = assemble(L_u_tent_mc)
        for bc_mc in u_bc_object_list_mc: bc_mc.apply(A_u_tent_mc, b_mc)
        solve(A_u_tent_mc, ut_mc.vector(), b_mc)

        # Compute pressure correction step.
        b_mc = assemble(L_p_corr_mc)
        for bc_mc in eta_bc_object_list_mc: bc_mc.apply(b_mc)
        if linear_divergence:
            a_p_corr_solver_mc.solve(eta1_mc.vector(), b_mc)
        else:
            A_p_corr_mc = assemble(a_p_corr_mc)
            for bc_mc in eta_bc_object_list_mc: bc_mc.apply(A_p_corr_mc)
            solve(A_p_corr_mc, eta1_mc.vector(), b_mc)

        # Compute Velocity correction step.
        b_mc = assemble(L_u_corr_mc)
        for bc_mc in u_bc_object_list_mc: bc_mc.apply(b_mc)
        a_u_corr_solver_mc.solve(u1_mc.vector(), b_mc)

        # Rotate functions for next time_step
        u00_mc.assign(u0_mc)
        u0_mc.assign(u1_mc)
        eta0_mc.assign(eta1_mc)

        # save to a result list
        result_u1.append(u1_mc(testnodes[0][0], testnodes[0][1])[0])
        result_v1.append(u1_mc(testnodes[0][0], testnodes[0][1])[1])
        result_eta1.append(eta1_mc(testnodes[0][0], testnodes[0][1]))


    return result_u1, result_v1, result_eta1


binrandom_u1_true = np.zeros([nsample, nsample, ntimestep+1])
binrandom_v1_true = np.zeros([nsample, nsample, ntimestep+1])
binrandom_eta1_true = np.zeros([nsample, nsample, ntimestep+1])

sample_initial_eta_key = initial_eta.keys()
sample_initial_u_key = initial_u.keys()
sample_boundary_u_key = boundary_u.keys()
sample_boundary_eta_key = boundary_eta.keys()
sample_bathy_key = bathymetry.keys()

sample_initial_eta_value = initial_eta.values()
sample_initial_u_value = initial_u.values()
sample_boundary_u_value = boundary_u.values()
sample_boundary_eta_value = boundary_eta.values()
sample_bathy_value = bathymetry.values()

for j in range(nsample):
    for k in range(nsample):
        # recal sto_viscosity sto_force
        q0 = samplex[j]
        q1 = sampley[k]
        x = sp.Symbol('x')
        y = sp.Symbol('y')
        sample_viscosity = str(eval(sto_viscosity)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        sample_friction = str(eval(sto_friction)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        sample_force_x = str(eval(sto_force_x)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        sample_force_y = str(eval(sto_force_y)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        # if sample_bathy.has_key('vary'):
        #     temp = sample_bathy.get('vary')
        #     temp1 = str(eval(temp)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        #     sample_bathy['vary'] = temp1
        # elif sample_bathy.has_key('class') and sample_bathy.values()[0][0] == 'type1':
        #     temp = sample_bathy.get('class')
        #     temp1 = str(eval(temp[3])).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        #     temp2 = str(eval(temp[4])).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        #     sample_bathy['class'][3] = temp1
        #     sample_bathy['class'][4] = temp2
        # else:
        #     pass
        # if sample_initial_eta.has_key('vary'):
        #     temp = initial_eta.get('vary')
        #     temp1 = str(eval(temp)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        #     sample_initial_eta['vary'] = temp1
        # if sample_initial_u.has_key('vary'):
        #     temp = initial_u.get('vary')
        #     temp1 = str(eval(temp[0])).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        #     temp2 = str(eval(temp[1])).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        #     sample_initial_u['vary'] = (temp1, temp2)
        # for i, key in enumerate(sample_boundary_u.keys()):
        #     value = sample_boundary_u[key]
        #     if isinstance(value, tuple) and isinstance(value[0], str):
        #         temp1 = str(eval(value[0])).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        #         temp2 = str(eval(value[1])).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        #         sample_boundary_u[key][i] = (temp1, temp2)
        # for i, key in enumerate(sample_boundary_eta.keys()):
        #     value = sample_boundary_eta[key]
        #     if isinstance(value, str) and value != 'freeslipyy' and value != 'freeslipxx' and value != 'noslip':
        #         temp = str(eval(value)).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        #         sample_boundary_eta[key][i] = temp
        sample_bathy = dict(zip(sample_bathy_key, sample_bathy_value))
        sample_initial_eta = dict(zip(sample_initial_eta_key, sample_initial_eta_value))
        sample_initial_u = dict(zip(sample_initial_u_key, sample_initial_u_value))
        sample_boundary_eta = dict(zip(sample_boundary_eta_key, sample_boundary_eta_value))
        sample_boundary_u = dict(zip(sample_boundary_u_key, sample_boundary_u_value))
        
        # test2_enable:
        #sample_initial_eta_value = str(eval(initial_eta.get('vary'))).replace("cos", "sp.cos").replace("sin", "sp.sin").replace("tan", "sp.tan")
        #sample_initial_eta = {sample_initial_eta_key[0]: sample_initial_eta_value}
        
        # test3_enable:
        # sample_bathy_value = bathymetry.get('class')[0:3]
        # temp3 = str(eval(bathymetry.get('class')[3]))
        # temp4 = str(eval(bathymetry.get('class')[4]))
        # sample_bathy_value.append(temp3)
        # sample_bathy_value.append(temp4)
        # sample_bathy = {sample_bathy_key[0]: sample_bathy_value}

        # test4_enable:
        # none here, all above

        # give u1 and eta1 at each time step here.
        result_u1_list, result_v1_list, result_eta1_list = func_mc(sample_viscosity, sample_friction, sample_force_x, sample_force_y, sample_bathy, sample_initial_u, sample_initial_eta, sample_boundary_u, sample_boundary_eta)
        for i in range(ntimestep+1):
            binrandom_u1_true[j, k, i] = result_u1_list[i]
            binrandom_v1_true[j, k, i] = result_v1_list[i]
            binrandom_eta1_true[j, k, i] = result_eta1_list[i]

        print "MC finish: " + str(j) + str(k)

# save MC and PC to file here.
np.save(outputdir+"binrandom_eta1_order_"+str(sto_poly_deg), binrandom_eta1)
np.save(outputdir+"binrandom_eta1_true_order_"+str(sto_poly_deg), binrandom_eta1_true)
np.save(outputdir+"binrandom_u1_order_"+str(sto_poly_deg), binrandom_u1)
np.save(outputdir+"binrandom_u1_true_order_"+str(sto_poly_deg), binrandom_u1_true)
np.save(outputdir+"binrandom_v1_order_"+str(sto_poly_deg), binrandom_v1)
np.save(outputdir+"binrandom_v1_true_order_"+str(sto_poly_deg), binrandom_v1_true)

for i in range(ntimestep+1):

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(sampleX, sampleY, binrandom_eta1[:, :, i], linewidth=0, antialiased=False)
    #fig.colorbar(surf, shrink=0.5, aspect=10)
    ax.scatter(sampleX, sampleY, binrandom_eta1_true[:, :, i], cmap=cm.cool, s=scatterplotsize, linewidth=0, antialiased=False)
    ax.set_xlim(coefficient[0], coefficient[1])
    ax.set_ylim(coefficient[2], coefficient[3])
    ###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
    #ax.set_zlim(zmin, zmax)
    ax.set_xticks(np.linspace(coefficient[0], coefficient[1], xyticknumber))
    ax.set_yticks(np.linspace(coefficient[2], coefficient[3], xyticknumber))
    ###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
    #ax.set_zticks(np.linspace(zmin, zmax, xyticknumber))
    ax.set_xlabel(r'$\xi_1$')
    ax.set_ylabel(r'$\xi_2$')
    ax.set_zlabel('ele')
    plt.title("(" + str(testnodes[0][0]) + ", " + str(testnodes[0][1]) + "); " + "stochastic order=" + str(sto_poly_deg) + "; time=" + str(i*timestep) + "s")
    #plt.show()
    plt.savefig(outputdir+"PCMC_3Dsurface_eta_point_" + str(testnodes[0][0]) + "_" + str(testnodes[0][1]) + "_stoorder_" + str(sto_poly_deg) + "_" + "{:02d}".format(i) + ".png")
    plt.close()
    
    temp = binrandom_eta1[:, :, i]-binrandom_eta1_true[:, :, i]
    plt.figure()
    plt.pcolor(sampleX, sampleY, temp, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    cbar = plt.colorbar(format='%.2e')
    cbar.ax.set_title(r'$eta_{pc}-eta_{mc}$')
    plt.xlim(coefficient[0], coefficient[1])
    plt.ylim(coefficient[2], coefficient[3])
    plt.xticks(np.linspace(coefficient[0], coefficient[1], xyticknumber))
    plt.yticks(np.linspace(coefficient[2], coefficient[3], xyticknumber))
    plt.xlabel(r'$\xi_1$')
    plt.ylabel(r'$\xi_2$')
    plt.title("(" + str(testnodes[0][0]) + ", " + str(testnodes[0][1]) + "); " + "stochastic order=" + str(sto_poly_deg) + "; time=" + str(i*timestep) + "s")
    #plt.show()
    plt.savefig(outputdir+"PCMC_2Ddifference_eta_point_" + str(testnodes[0][0]) + "_" + str(testnodes[0][1]) + "_stoorder_" + str(sto_poly_deg) + "_" + "{:02d}".format(i) + ".png")
    plt.close()
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(sampleX, sampleY, binrandom_u1[:, :, i], linewidth=0, antialiased=False)
    #fig.colorbar(surf, shrink=0.5, aspect=10)
    ax.scatter(sampleX, sampleY, binrandom_u1_true[:, :, i], cmap=cm.cool, s=scatterplotsize, linewidth=0, antialiased=False)
    ax.set_xlim(coefficient[0], coefficient[1])
    ax.set_ylim(coefficient[2], coefficient[3])
    ###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
    #ax.set_zlim(zmin, zmax)
    ax.set_xticks(np.linspace(coefficient[0], coefficient[1], xyticknumber))
    ax.set_yticks(np.linspace(coefficient[2], coefficient[3], xyticknumber))
    ###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
    #ax.set_zticks(np.linspace(zmin, zmax, xyticknumber))
    ax.set_xlabel(r'$\xi_1$')
    ax.set_ylabel(r'$\xi_2$')
    ax.set_zlabel('u')
    plt.title("(" + str(testnodes[0][0]) + ", " + str(testnodes[0][1]) + "); " + "stochastic order=" + str(sto_poly_deg) + "; time=" + str(i*timestep) + "s")
    #plt.show()
    plt.savefig(outputdir+"PCMC_3Dsurface_u_point_" + str(testnodes[0][0]) + "_" + str(testnodes[0][1]) + "_stoorder_" + str(sto_poly_deg) + "_" + "{:02d}".format(i) + ".png")
    plt.close()
    
    temp = binrandom_u1[:, :, i]-binrandom_u1_true[:, :, i]
    plt.figure()
    plt.pcolor(sampleX, sampleY, temp, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    cbar = plt.colorbar(format='%.2e')
    cbar.ax.set_title(r'$u_{pc}-u_{mc}$')
    plt.xlim(coefficient[0], coefficient[1])
    plt.ylim(coefficient[2], coefficient[3])
    plt.xticks(np.linspace(coefficient[0], coefficient[1], xyticknumber))
    plt.yticks(np.linspace(coefficient[2], coefficient[3], xyticknumber))
    plt.xlabel(r'$\xi_1$')
    plt.ylabel(r'$\xi_2$')
    plt.title("(" + str(testnodes[0][0]) + ", " + str(testnodes[0][1]) + "); " + "stochastic order=" + str(sto_poly_deg) + "; time=" + str(i*timestep) + "s")
    #plt.show()
    plt.savefig(outputdir+"PCMC_2Ddifference_u_point_" + str(testnodes[0][0]) + "_" + str(testnodes[0][1]) + "_stoorder_" + str(sto_poly_deg) + "_" + "{:02d}".format(i) + ".png")
    plt.close()
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    surf = ax.plot_surface(sampleX, sampleY, binrandom_v1[:, :, i], linewidth=0, antialiased=False)
    #fig.colorbar(surf, shrink=0.5, aspect=10)
    ax.scatter(sampleX, sampleY, binrandom_v1_true[:, :, i], cmap=cm.cool, s=scatterplotsize, linewidth=0, antialiased=False)
    ax.set_xlim(coefficient[0], coefficient[1])
    ax.set_ylim(coefficient[2], coefficient[3])
    ###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
    #ax.set_zlim(zmin, zmax)
    ax.set_xticks(np.linspace(coefficient[0], coefficient[1], xyticknumber))
    ax.set_yticks(np.linspace(coefficient[2], coefficient[3], xyticknumber))
    ###### WARNING!!! CHANGE IT if you choose the same sto_order but have MORE THAN ONE NODEs.
    #ax.set_zticks(np.linspace(zmin, zmax, xyticknumber))
    ax.set_xlabel(r'$\xi_1$')
    ax.set_ylabel(r'$\xi_2$')
    ax.set_zlabel('v')
    plt.title("(" + str(testnodes[0][0]) + ", " + str(testnodes[0][1]) + "); " + "stochastic order=" + str(sto_poly_deg) + "; time=" + str(i*timestep) + "s")
    #plt.show()
    plt.savefig(outputdir+"PCMC_3Dsurface_v_point_" + str(testnodes[0][0]) + "_" + str(testnodes[0][1]) + "_stoorder_" + str(sto_poly_deg) + "_" + "{:02d}".format(i) + ".png")
    plt.close()
    
    temp = binrandom_v1[:, :, i]-binrandom_v1_true[:, :, i]
    plt.figure()
    plt.pcolor(sampleX, sampleY, temp, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    cbar = plt.colorbar(format='%.2e')
    cbar.ax.set_title(r'$v_{pc}-v_{mc}$')
    plt.xlim(coefficient[0], coefficient[1])
    plt.ylim(coefficient[2], coefficient[3])
    plt.xticks(np.linspace(coefficient[0], coefficient[1], xyticknumber))
    plt.yticks(np.linspace(coefficient[2], coefficient[3], xyticknumber))
    plt.xlabel(r'$\xi_1$')
    plt.ylabel(r'$\xi_2$')
    plt.title("(" + str(testnodes[0][0]) + ", " + str(testnodes[0][1]) + "); " + "stochastic order=" + str(sto_poly_deg) + "; time=" + str(i*timestep) + "s")
    #plt.show()
    plt.savefig(outputdir+"PCMC_2Ddifference_v_point_" + str(testnodes[0][0]) + "_" + str(testnodes[0][1]) + "_stoorder_" + str(sto_poly_deg) + "_" + "{:02d}".format(i) + ".png")
    plt.close()
