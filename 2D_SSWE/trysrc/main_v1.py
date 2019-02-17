from fenics import *
import sys
from input import *
from makestobasis import MakeStoBasis
from makestomodes import MakeStoModes
from makestoijk import MakeStoIJK
from makemesh import MakeMesh
from make_bath import make_bath
from makeic import MakeInitialObjectFunc
from makebc import MakeBoundaryObjectList

set_log_level(INFO)
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
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

# Get equation parameter settingTrue
g = Constant(gravity)

# Get the bathymetry
project_bottom_function = make_bath(bathymetry, basis_struct, P)

# Get initial condition in stochastic setting
u_ic_function, eta_ic_function = MakeInitialObjectFunc(initial_u, initial_eta, basis_struct, W, P)

# Get boundary condition in stochastic setting
u_bc_object_list, eta_bc_object_list, uTimeDependent, etaTimeDependent, u_list_expression, eta_list_expression = MakeBoundaryObjectList(boundary_u, boundary_eta, basis_struct, bcfile, mesh, domain, W, P, t)



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

A = FunctionSpace(mesh, V)
B = FunctionSpace(mesh, Q)
f_u_expression_object_list = []
f_v_expression_object_list = []
friction_expression_object_list = []
nu_expression_object_list = []
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
    nu_expression_object_list.append(Expression(nu_expression_list[k], element=B.ufl_element()))

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
                F_u_tent = F_u_tent + ( stoIJK[i][j][k] * friction_expression_object_list[i] * inner(u_mean[2*j], v[2*k]) * dx +
                                        stoIJK[i][j][k] * friction_expression_object_list[i] * inner(u_mean[2*j+1], v[2*k+1]) * dx )
                if include_viscosity:
                    F_u_tent = F_u_tent + ( stoIJK[i][j][k] * nu_expression_object_list[i] * inner(grad(v[2*k]), grad(u_mean[2*j])) * dx +
                                            stoIJK[i][j][k] * nu_expression_object_list[i] * inner(grad(v[2*k+1]), grad(u_mean[2*j+1])) * dx )
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
                    F_p_corr = F_p_corr + ( stoIJK[i][j][k] * theta**2 * g * dt**2 * H * inner(grad(q), grad(eta_diff)) * dx
                                    + stoIJK[i][j][k] * dt * q * grad(H*ut_mean[0])[0] * dx
                                    + stoIJK[i][j][k] * dt * q * grad(H*ut_mean[1])[1] * dx )
                else:
                    F_p_corr = F_p_corr + ( stoIJK[i][j][k] * theta**2 * g * dt**2 * H[i] * inner(grad(q[k]), grad(eta_diff[j])) * dx
                                    + stoIJK[i][j][k] * dt * q[k] * grad(H[i]*ut_mean[2*j])[0] * dx
                                    + stoIJK[i][j][k] * dt * q[k] * grad(H[i]*ut_mean[2*j+1])[1] * dx )
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

if nmodes == 1:
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
        u_file.append(XDMFFile("u_"+outputstr+"{:02d}".format(mode)+".xdmf"))
        eta_file.append(XDMFFile("eta_"+outputstr+"{:02d}".format(mode)+".xdmf"))
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

