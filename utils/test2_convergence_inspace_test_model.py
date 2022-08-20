" Write a nonlinear shallow water equation solver."
from dolfin.cpp.log import log
from fenics import *
from mshr import *
import numpy as np
import sympy as sp



hsize = 16
# loc = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/test2_res/20181001_generalized/mode_0_fine_mesh/convergence_test_result/10mesh_5mesh_small_time_step/"
loc = "./"

x0 = 0
y0 = 0
x1 = 100
y1 = 50
nx = 10*hsize
ny = 5*hsize
startime = 0.0
endtime = 0.01
timestep = 0.01
theta = 0.5
gravity = 9.81
depth = 20.0
viscosity = 0.0
friction = 0.0

# error norm in the whole space integral

aa = 0.1
LL = 100
gg = 9.81
HH = 20
pi = 3.1416



# build up mesh
mesh = RectangleMesh(Point(x0, y0), Point(x1, y1), nx, ny)


# Get Funciton space
V = VectorFunctionSpace(mesh, 'CG', 2, dim=2)
Q = FunctionSpace(mesh, 'CG', 1)


# Get temporal setting
theta = Constant(theta)
dt = Constant(timestep)
finish_time = Constant(endtime)
t = Constant(startime)


# Get equation setting
g = Constant(gravity)
h = Constant(depth)
nu = Constant(viscosity)
friction = Constant(friction)
include_viscosity = True
include_advection = True
linear_divergence = False
f_u = Constant((0, 0))
initial_condition_u = Constant((0, 0))
initial_condition_eta = Expression(" 0.1 * cos(pi/100.0*x[0]) ", element=Q.ufl_element())
# initial_condition_eta = Constant(0.0)

# Get boundary condition ==> no neumann boundary condition
def leftbc(x,on_boundary):
    return on_boundary and near(x[0], x0)
def rightbc(x,on_boundary):
    return on_boundary and near(x[0], x1)
def updownbc(x,on_boundary):
    return on_boundary and (near(x[1], y0) or near(x[1], y1))

# Get boundary normal direction
n = FacetNormal(mesh)

u_left =  Constant(0)
u_right = Constant(0)
u_updown = Constant(0)
bc_left = DirichletBC(V.sub(0), u_left, leftbc, method="topological")
bc_right = DirichletBC(V.sub(0), u_right, rightbc, method="topological")
bc_updown = DirichletBC(V.sub(1), u_updown, updownbc, method="topological")

bcu = [bc_left, bc_right, bc_updown]
bceta = []

# Get Test and Trial Functions
v = TestFunction(V)
u = TrialFunction(V)
q = TestFunction(Q)
eta = TrialFunction(Q)


# Functions
u00 = Function(V)
u0 = Function(V, name="u0")
ut = Function(V) # Tentative velocity
u1 = Function(V, name="u")
eta0 = Function(Q, name="eta0")
eta1 = Function(Q, name="eta")


# Define the water depth
if linear_divergence:
    H = h
else:
    H = eta0 + h


# Load initial condition
u_ic = project(initial_condition_u, V)
u0.assign(u_ic)
u00.assign(u_ic)
u1.assign(u_ic)

eta_ic = project(initial_condition_eta, Q)
eta0.assign(eta_ic)
eta1.assign(eta_ic)


# Tentative Velocity step
u_mean = theta * u + (1. - theta) * u0
u_bash = 3./2 * u0 - 1./2 * u00
u_diff = u - u0
norm_u0 = inner(u0, u0)**0.5
F_u_tent = ((1/dt) * inner(v, u_diff) * dx +
            g * inner(v, grad(eta0)) * dx +
            friction / H * norm_u0 * inner(u_mean, v) * dx -
            inner(v, f_u) * dx )

if include_viscosity:
    F_u_tent += nu * inner(grad(v), grad(u_mean)) * dx

if include_advection:
    F_u_tent += inner(v, grad(u_mean)*u_bash) * dx
a_u_tent = lhs(F_u_tent)
L_u_tent = rhs(F_u_tent)


# Pressure correction step
eta_diff = eta - eta0
ut_mean = theta * ut + (1. - theta) * u0
F_p_corr = (q*eta_diff + g * dt**2 * theta**2 * H * inner(grad(q), grad(eta_diff)))*dx + dt*q*div(H*ut_mean)*dx

a_p_corr = lhs(F_p_corr)
L_p_corr = rhs(F_p_corr)


# Velocity correction step
eta_diff = eta1 - eta0
a_u_corr = inner(v, u)*dx
L_u_corr = inner(v, ut)*dx - dt*g*theta*inner(v, grad(eta_diff))*dx


# Assemble matrices
A_u_corr = assemble(a_u_corr)
for bc in bcu: bc.apply(A_u_corr)
a_u_corr_solver = LUSolver(A_u_corr)
a_u_corr_solver.parameters["reuse_factorization"] = True


if linear_divergence:
    A_p_corr = assemble(a_p_corr)
    for bc in bceta: bc.apply(A_p_corr)
    a_p_corr_solver = LUSolver(A_p_corr)
    a_p_corr_solver.parameters["reuse_factorization"] = True


# Time stepping.
log(INFO, "Start of time loop")
timestep = 0
u_file = XDMFFile(loc+"u_"+str(hsize)+".xdmf")
eta_file = XDMFFile(loc+"eta_"+str(hsize)+".xdmf")
u_file.write(u1, float(t))
eta_file.write(eta1, float(t))
error_u_l2 = []
error_eta_l2 = []
error_u_h1 = []
error_eta_h1 = []


while float(t - finish_time) <= - 1e3*DOLFIN_EPS:

    # Update time_step
    timestep += 1
    t = Constant(t+dt)


    # Update bc's
    t_theta = Constant(t - (1.0 - theta) * dt)


    # Compute tentative velocity step
    log(INFO, "Solve for tentative velocity.")
    A_u_tent = assemble(a_u_tent)
    b = assemble(L_u_tent)
    for bc in bcu: bc.apply(A_u_tent, b)

    solve(A_u_tent, ut.vector(), b)


    # Compute pressure correction step
    log(INFO, "Solve for pressure correction.")
    b = assemble(L_p_corr)
    for bc in bceta: bc.apply(b)

    if linear_divergence:
        a_p_corr_solver.solve(eta1.vector(), b)
    else:
        A_p_corr = assemble(a_p_corr)
        for bc in bceta: bc.apply(A_p_corr)
        solve(A_p_corr, eta1.vector(), b)


    # Compute Velocity correction step
    log(INFO, "Solve for velocity update.")
    b = assemble(L_u_corr)
    for bc in bcu: bc.apply(b)

    a_u_corr_solver.solve(u1.vector(), b)

    # Write to file
    u_file.write(u1, float(t))
    eta_file.write(eta1, float(t))

    # Rotate functions for next time_step
    u00.assign(u0)
    u0.assign(u1)
    eta0.assign(eta1)

    temp1 = Expression(["aa * pow(gg*HH, 0.5) / HH * sin(pi / LL * x[0]) * sin(pi * pow(gg*HH, 0.5) / LL * t)", "0"], aa=aa, gg=gg, HH=HH, LL=LL, pi=pi, t=float(t), degree=5)
    temp2 = Expression("aa * cos(pi / LL * x[0]) * cos(pi * pow(gg*HH, 0.5) / LL * t)", aa=aa, gg=gg, HH=HH, LL=LL, pi=pi, t=float(t), degree=4)

    trueufuc = project(temp1, V)
    trueetafuc = project(temp2, Q)

    error_u_l2.append(errornorm(trueufuc, u1, "L2", 5))
    error_eta_l2.append(errornorm(trueetafuc, eta1, "L2", 4))

    error_u_h1.append(errornorm(trueufuc, u1, "H1", 5))
    error_eta_h1.append(errornorm(trueetafuc, eta1, "H1", 4))

    print(t)

# f = open(loc+'error_u_l2_'+str(hsize), "w+")
# for err in error_u_l2:
#     f.write(str(err))
#     f.write("\n")
# f.close()

f2 = open(loc+'error_eta_l2_'+str(hsize), "w+")
for err in error_eta_l2:
    f2.write(str(err))
    f2.write("\n")
f2.close()

# f = open(loc+'error_u_h1_'+str(hsize), "w+")
# for err in error_u_h1:
#     f.write(str(err))
#     f.write("\n")
# f.close()

# f2 = open(loc+'error_eta_h1_'+str(hsize), "w+")
# for err in error_eta_h1:
#     f2.write(str(err))
#     f2.write("\n")
# f2.close()
