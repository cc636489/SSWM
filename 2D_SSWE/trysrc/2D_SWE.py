" Write a nonlinear shallow water equation solver."


from fenics import *
from mshr import *
import numpy as np
import sympy as sp
import pdb
import matplotlib.pyplot as plt


set_log_level(INFO)
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}

startime = 0.0
endtime = 600.0
timestep = 5.0
theta = 0.5
gravity = 9.81
viscosity = 30
friction = 0.003
include_viscosity = True
include_advection = True
linear_divergence = False
include_les = True
les_parameters = {'smagorinsky_coefficient': 1e-2}


mesh = Mesh("inlet.xml")
V = VectorFunctionSpace(mesh, 'CG', 2, dim=2)
Q = FunctionSpace(mesh, 'CG', 1)
v = TestFunction(V)
u = TrialFunction(V)
q = TestFunction(Q)
eta = TrialFunction(Q)
u00 = Function(V)
u0 = Function(V, name="u0")
ut = Function(V) # Tentative velocity
u1 = Function(V, name="u")
eta0 = Function(Q, name="eta0")
eta1 = Function(Q, name="eta")



# Get temporal setting
theta = Constant(theta)
dt = Constant(timestep)
finish_time = Constant(endtime)
t = Constant(startime)

# Get equation constant
g = Constant(gravity)
nu = Constant(viscosity)
friction = Constant(friction)

# Get source term setting
f_u = Constant((0.0, 0))

# Get bottom shape
class bottomExpression(Expression):
    def eval(self, value, x):
        if x[0] <= 2250:
            value[0] = - 14.0/2250.0 * x[0] + 19.0
        else:
            value[0] = 5.0
    def value_shape(self):
        return (1,)
bottom = bottomExpression(element=Q.ufl_element())
project_bottom = interpolate(bottom, Q)
if linear_divergence:
    H = project_bottom
else:
    H = eta0 + project_bottom

# Get initial condition
initial_condition_u = Constant((0, 0))
initial_condition_eta = Constant(0)
u_ic = interpolate(initial_condition_u, V)
u0.assign(u_ic)
u00.assign(u_ic)
u1.assign(u_ic)
eta_ic = interpolate(initial_condition_eta, Q)
eta0.assign(eta_ic)
eta1.assign(eta_ic)

# Get boundary condition
boundaries = MeshFunction("size_t", mesh, "inlet_facet_region.xml")
#opentidal = Expression(("sin(pi*t/5)*x[1]*(3000.0-x[1])/2250000.0", "0.0"), element = V.ufl_element(), t=t)
opentidal = Expression(" 0.5 * sin(pi/100.0 * t)", element = Q.ufl_element(), t=t)
bc_open = DirichletBC(Q, opentidal, boundaries, 3)
bc_noslip_xx = DirichletBC(V.sub(0), Constant(0), boundaries, 1)
bc_noslip_yy = DirichletBC(V.sub(1), Constant(0), boundaries, 2)
bcu = [bc_noslip_xx, bc_noslip_yy]
bceta = [bc_open]

# For Large eddy model
class LES(object):

    def __init__(self, V, u, smagorinsky_coefficient):
        self._V = V
        self.eddy_viscosity = Function(V)

        # Create a eddy viscosity solver
        les_lhs, les_rhs = self._eddy_viscosity_eqn(u, smagorinsky_coefficient)

        eddy_viscosity_problem = LinearVariationalProblem(les_lhs, les_rhs,
                 self.eddy_viscosity, bcs=[])
        self._solver = LinearVariationalSolver(eddy_viscosity_problem)
        self._solver.parameters["linear_solver"] = "lu"
        self._solver.parameters["symmetric"] = True
        self._solver.parameters["lu_solver"]["reuse_factorization"] = True

    def _strain_rate_tensor(self, u):
        S = 0.5*(grad(u) + grad(u).T)
        return S

    def _eddy_viscosity_eqn(self, u, smagorinsky_coefficient):

        dim = len(u)
        w = TestFunction(self._V)
        eddy_viscosity = TrialFunction(self._V)

        cell_vol = CellVolume(self._V.mesh())
        filter_width = cell_vol**(1.0/dim)

        S = self._strain_rate_tensor(u)
        second_invariant = 0.0
        for i in range(0, dim):
           for j in range(0, dim):
              second_invariant += 2.0*(S[i,j]**2)

        second_invariant = sqrt(second_invariant)
        rhs = (smagorinsky_coefficient*filter_width)**2*second_invariant

        lhs = inner(w, eddy_viscosity)*dx
        rhs = inner(w, rhs)*dx

        return lhs, rhs

    def solve(self):
        self._solver.solve()
        return self.eddy_viscosity
if include_les:
	les_V = FunctionSpace(mesh, "CG", 1)
	les = LES(les_V, u0, les_parameters['smagorinsky_coefficient'])
	eddy_viscosity = les.eddy_viscosity
	nu = Constant(viscosity)+ eddy_viscosity
	project_nu = project(nu, les_V)
else:
	eddy_viscosity = None

#######################################################################
#######  Solving the equation starts ...
#######################################################################


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
    F_u_tent += inner(v, grad(u_mean)*u_bash) * dx ## FIXME!! need concern about the actual splitting form. should be grad(u_mean)*u_bash.
a_u_tent = lhs(F_u_tent)
L_u_tent = rhs(F_u_tent)


# Pressure correction step
eta_diff = eta - eta0
ut_mean = theta * ut + (1. - theta) * u0
F_p_corr = (q*eta_diff + g * dt**2 * theta**2 * H * inner(grad(q), grad(eta_diff)))*dx + dt*q*div( H * ut_mean)*dx
#pdb.set_trace()
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
u_file = XDMFFile("u.xdmf")
eta_file = XDMFFile("eta.xdmf")
eddy_vis = XDMFFile("nu.xdmf")
u_file.write(u1, float(t))
eta_file.write(eta1, float(t))
eddy_vis.write(project_nu, float(t))

while float(t - finish_time) <= - 1e3*DOLFIN_EPS:

	timestep += 1
	t = Constant(t+dt)
	print timestep

	# update source term
	t_theta = Constant(t - (1.0 - theta) * dt)

	# update boundary condition
	opentidal.user_parameters['t'] = t

	# update eddy_viscosity
	if include_les:
		log(INFO, "Compute eddy viscosity.")
		les.solve()
		eddy_viscosity = les.eddy_viscosity
		nu = Constant(viscosity) + eddy_viscosity
		project_nu = project(nu, les_V)
#F_u_tent=((1/dt)*inner(v,u_diff)*dx+g*inner(v,grad(eta0))*dx+friction/H*u0_norm*inner(u_mean,v)*dx-inner(v,f_u)*dx)+nu*inner(grad(v), grad(u_mean)) * dx

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
	eddy_vis.write(project_nu, float(t))

#	plot(u1, interactive=True)
#	plot(eta1, interactive=True)

	# Rotate functions for next time_step
	u00.assign(u0)
	u0.assign(u1)
	eta0.assign(eta1)


log(INFO, "End of time loop.")
print "time_step is ", timestep