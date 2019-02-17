" Write a nonlinear shallow water equation solver."


from fenics import *
from mshr import *
import numpy as np
import sympy as sp
import pdb


set_log_level(INFO)
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}

hsize = [2, 4, 8, 16]
loc = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/test3_res_wiLES/20181001_generalized/convergence_test/Initial_time_step_0.01_second/'


def calculation(hsize):
    x0 = 0
    y0 = 0
    x1 = 1000
    y1 = 200
    nx = 10*hsize
    ny = 2*hsize
    startime = 0.0
    endtime = 0.01
    timestep = 0.01
    theta = 0.5
    gravity = 9.81
    viscosity = 0.0
    friction = 0.0
    
    
    # Get equation constant
    g = Constant(gravity)
    nu = Constant(viscosity)
    friction = Constant(friction)
    include_viscosity = True
    include_advection = True
    linear_divergence = False
    include_les = False
    les_parameters = {'smagorinsky_coefficient': 0.0}
    
    
    # Get force setting
    f_u = Constant((0, 0))
    
    
    # Get temporal setting
    theta = Constant(theta)
    dt = Constant(timestep)
    finish_time = Constant(endtime)
    t = Constant(startime)
    
    
    # Get initial condition
    initial_condition_u = Constant((0, 0))
    initial_condition_eta = Constant(0)
    
    
    # build up mesh
    mesh = RectangleMesh(Point(x0, y0), Point(x1, y1), nx, ny)
    
    
    # Get Funciton space
    V = VectorFunctionSpace(mesh, 'CG', 2, dim=2)
    Q = FunctionSpace(mesh, 'CG', 1)
    
    
    # Get Test and Trial Functions
    v = TestFunction(V)
    u = TrialFunction(V)
    q = TestFunction(Q)
    eta = TrialFunction(Q)
    
    
    
    # Get boundary condition ==> no neumann boundary condition
    def leftbc(x,on_boundary):
        return on_boundary and near(x[0], x0)
    def rightbc(x,on_boundary):
        return on_boundary and near(x[0], x1)
    def updownbc(x,on_boundary):
        return on_boundary and (near(x[1], y0) or near(x[1], y1))
    
    u_left = Constant(0)
    eta_right = Expression(" 0.1 * sin(pi/20.0 * t) ", element=Q.ufl_element(), t=t)
    u_updown = Constant(0)
    bc_left = DirichletBC(V.sub(0), u_left, leftbc, method="topological")
    bc_right = DirichletBC(Q, eta_right, rightbc)
    bc_updown = DirichletBC(V.sub(1), u_updown, updownbc, method="topological")
    
    bcu = [bc_left, bc_updown]
    bceta = [bc_right]
    
    
    
    # Functions
    u00 = Function(V)
    u0 = Function(V, name="u0")
    ut = Function(V) # Tentative velocity
    u1 = Function(V, name="u")
    eta0 = Function(Q, name="eta0")
    eta1 = Function(Q, name="eta")
    
    
    # Load initial condition
    u_ic = interpolate(initial_condition_u, V)
    u0.assign(u_ic)
    u00.assign(u_ic)
    u1.assign(u_ic)
    
    eta_ic = interpolate(initial_condition_eta, Q)
    eta0.assign(eta_ic)
    eta1.assign(eta_ic)
    
    
    # Get bottom shape
    class bottomExpression(Expression):
    
        def eval(self, value, x):
            if x[0] < 400 or x[0] > 600:
                value[0] = 5.0
            else:
                value[0] = -3*((x[0]-500)/100)**4 + 6*((x[0]-500)/100)**2 + 2
    
        def value_shape(self):
            return (1, )
    
    # Define the water depth
    bottom = bottomExpression(element=Q.ufl_element())
    # bottom = Expression("20.0+x[0]+x[1]", degree=1, cell=triangle)
    # bottom =  Constant(20)
    project_bottom = interpolate(bottom, Q)
    if linear_divergence:
        H = project_bottom
    else:
        H = eta0 + project_bottom
    
    
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
        nu = Constant(viscosity) + eddy_viscosity
        project_nu = project(nu, les_V)
    else:
        eddy_viscosity = None
    
    
    
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
    F_p_corr = (q*eta_diff + g * dt**2 * theta**2 * H * inner(grad(q), grad(eta_diff)))*dx + dt*q*div( H * ut_mean)*dx
    #pdb.set_trace()
    a_p_corr = lhs(F_p_corr)
    L_p_corr = rhs(F_p_corr)
    
    
    # Velocity correction step
    eta_diff = eta1 - eta0
    F_u_corr = inner(v, u)*dx - inner(v, ut)*dx + dt*g*theta*inner(v, grad(eta_diff))*dx
    a_u_corr = lhs(F_u_corr)
    L_u_corr = rhs(F_u_corr)
    #a_u_corr = inner(v, u)*dx
    #L_u_corr = inner(v, ut)*dx - dt*g*theta*inner(v, grad(eta_diff))*dx
    
    
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
    u_file = XDMFFile(loc+"u_woles.pvd")
    eta_file = XDMFFile(loc+"eta_woles.pvd")
    u_file.write(u1, float(t))
    eta_file.write(eta1, float(t))
    
    while float(t - finish_time) <= - 1e3*DOLFIN_EPS:
    
        # Update time_step
        timestep += 1
        t = Constant(t+dt)
        print timestep
        print "time is: ", float(t)
    
    
        # Update bc's
        t_theta = Constant(t - (1.0 - theta) * dt)
        eta_right.user_parameters['t'] = t
        bc_right = DirichletBC(Q, eta_right, rightbc)
    
        # update eddy_viscosity
        if include_les:
            log(INFO, "Compute eddy viscosity.")
            les.solve()
            eddy_viscosity = les.eddy_viscosity
            nu = Constant(viscosity) + eddy_viscosity
            project_nu = project(nu, les_V)
    
            # for i in range(0, 1010, 10):
            #     for j in range(0, 210, 10):
            #         print("i = ", i, "j = ", j, project_nu(i, j))
    
    
        # Compute tentative velocity step
        log(INFO, "Solve for tentative velocity.")
        A_u_tent = assemble(a_u_tent)
        # print("utent matrix:")
        # print(A_u_tent.getrow(0))
        b = assemble(L_u_tent)
        # print("utent vector:")
        # print(b.array())
        for bc in bcu: bc.apply(A_u_tent, b)
    
        solve(A_u_tent, ut.vector(), b)
    
    
        # Compute pressure correction step
        log(INFO, "Solve for pressure correction.")
        b = assemble(L_p_corr)
        # print("pcorr vector:")
        # print(b.array())
        for bc in bceta: bc.apply(b)
    
        if linear_divergence:
            a_p_corr_solver.solve(eta1.vector(), b)
        else:
            A_p_corr = assemble(a_p_corr)
            # print("pcorr matrix:")
            # print(A_p_corr.getrow(0))
            for bc in bceta: bc.apply(A_p_corr)
            solve(A_p_corr, eta1.vector(), b)
    
    
        # Compute Velocity correction step
        log(INFO, "Solve for velocity update.")
        b = assemble(L_u_corr)
        # print("ucorr vector:")
        # print(b.array()[-10:-1])
        for bc in bcu: bc.apply(b)
    
        # print("ucorr matrix:")
        # print(A_u_corr.getrow(0))
        a_u_corr_solver.solve(u1.vector(), b)
    
    
        # Write to file
        u_file.write(u1, float(t))
        eta_file.write(eta1, float(t))
        
    
    
        # Rotate functions for next time_step
        u00.assign(u0)
        u0.assign(u1)
        eta0.assign(eta1)
    
        # for i in range(500, 1010, 10):
        #     for j in range(0, 210, 10):
        #         print("i = ", i,"j = ", j, "u1 = ", u1(i, j))
        # for i in range(500, 1010, 10):
        #     for j in range(0, 210, 10):
        #         print("i = ", i,"j = ", j, "eta1 = ", eta1(i, j))
    
        
        print("===========================================")

        return u1, eta1, V, Q

u1, eta1, V, Q = calculation(hsize=64)
u_true, eta_true = Function(V), Function(Q)
u_true.vector()[:] = u1.vector().array()
eta_true.vector()[:] = eta1.vector().array()

error_u = []
error_eta = []
for h in hsize:
    u_model, eta_model, _, _ = calculation(hsize=h)
    error_u.append(errornorm(u_true, u_model, 'L2', 6))
    error_eta.append(errornorm(eta_true, eta_model, 'L2', 4))
print error_u
print error_eta

f = open(loc+'error_u_timestep1_woles', "w+")
for err in error_u:
    f.write(str(err))
    f.write("\n")
f.close()

f2 = open(loc+'error_eta_timestep1_woles', "w+")
for err in error_eta:
    f2.write(str(err))
    f2.write("\n")
f2.close()

