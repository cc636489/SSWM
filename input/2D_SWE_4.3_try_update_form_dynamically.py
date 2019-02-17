" Write a nonlinear shallow water equation solver."


from fenics import *
from netCDF4 import *
import sys

set_log_level(INFO)
parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, \
               "eliminate_zeros": True, \
               "precompute_basis_const": True, \
               "precompute_ip_const": True}
sys.setrecursionlimit(10000)

tidal_amplitude = 0.75  # galveston bay area normal tides.
tidal_period = 12.42*60*60
startime = 0.0
endtime = 223500
dt = 447
theta = 1.0
gravity = 9.81
viscosity = 1e-6
friction = 0.0
include_viscosity = True
include_advection = True
linear_divergence = False
include_les = False
include_supg = True
include_crosswind = True
include_ip = False
include_viscosity_dgu = False
les_parameters = {'smagorinsky_coefficient': 0.17}
prefix = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/" \
         "test8_SUPG_SamgLily_crosswind/test4_supg_crosswind_try_update_dynamic_woles_wofric"

mesh = Mesh("inlet_dgswem_compare_finer.xml")
V = VectorFunctionSpace(mesh, 'CG', 2, dim=2)
Q = FunctionSpace(mesh, 'CG', 1)
v = TestFunction(V)
u = TrialFunction(V)
q = TestFunction(Q)
eta = TrialFunction(Q)
u00 = Function(V)
u0 = Function(V, name="u0")
ut = Function(V)
u1 = Function(V, name="u")
eta0 = Function(Q, name="eta0")
eta1 = Function(Q, name="eta")


# Get temporal setting
theta = Constant(theta)
dt = Constant(dt)
finish_time = Constant(endtime)
t = Constant(startime)


# Get equation constant
g = Constant(gravity)
nu = Constant(viscosity)
friction = Constant(friction)


# Get source term setting
f_u = Constant((0.0, 0.0))


# Get bottom shape
code = '''

namespace dolfin {

struct node {
    double x;
    double y;
    double v;
    node *next;
};

class MyFun : public Expression
{
    private:
        node *root;

    public:

        MyFun(): Expression()
        {
            root = new node;
            root->next = 0;
        };

        void update(double _x, double _y, double _v)
        {
            node *newNode = new node;
            newNode->x = _x;
            newNode->y = _y;
            newNode->v = _v;
            newNode->next = root;
            root = newNode;
        };

        void eval(Array<double>& values, const Array<double>& x) const
        {
            double vv = findval(x[0], x[1]);
            values[0] = vv;
            //cout << x[0] << "  " << x[1] << "  " << vv << endl;
            //cout << "  " << endl;
        };

        double findval(double _x, double _y) const
        {
            node *p = root;
            while (p->next != 0)
            {
            // Assume that root node has biggest x-value.
            // Traverse down the list until p->x = _x and p->y = _y, then assign p->v to _v, and return value, break.
                if ( (fabs(p->x - _x)<1.0e-4) && (fabs(p->y - _y)<1.0e-4) )
                {
                    //cout << fabs(p->x-_x) << "  " << fabs(p->y-_y) << " " << p->x << " " << _x << " " << p->y << " " << _y << endl;
                    double find = p->v;
                    return find;
                }
                else
                {
                    p = p->next;
                }
            }
            return 0;

        };

};
};

'''
nc = Dataset("inlet_dgswem_compare_finer.nc", 'r', format='NETCDF4')
y_coord = nc.variables['y_coord'][:].data.tolist()
x_coord = nc.variables['x_coord'][:].data.tolist()
bath_depth = nc.variables['bathy'][:].data.tolist()
bottom = Expression(code, element=Q.ufl_element())
for j in range(len(x_coord)): bottom.update(x_coord[j], y_coord[j], bath_depth[j])
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
boundaries = MeshFunction("size_t", mesh, "inlet_dgswem_compare_finer_facet_region.xml")
opentidal = Expression("amp*sin(omega*t)", element=Q.ufl_element(), t=Constant(0), amp=tidal_amplitude, omega=2*pi/tidal_period, g=gravity, H=project_bottom)
bc_open = DirichletBC(Q, opentidal, boundaries, 2)
bcu = []
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
    les = LES(Q, u0, les_parameters['smagorinsky_coefficient'])
    eddy_viscosity = les.eddy_viscosity
    nu = Constant(viscosity) + eddy_viscosity
project_nu = project(nu, Q)

#######################################################################
#######  Solving the equation starts ...
#######################################################################

h = CellSize(mesh)
n = FacetNormal(mesh)

# Tentative Velocity step
u_mean = theta * u + (1. - theta) * u0
u0_mean = theta * u0 + (1. - theta) * u00
u_bash = 3./2 * u0 - 1./2 * u00
u_diff = u - u0
u_norm = inner(u0, u0) ** 0.5
du_norm = inner(grad(u0[0]), grad(u0[0])) ** 0.5
dv_norm = inner(grad(u0[1]), grad(u0[1])) ** 0.5

F_u_tent = ((1/dt) * inner(v, u_diff) * dx +
            g * inner(v, grad(eta0)) * dx +
            friction / H * u_norm * inner(u_mean, v) * dx -
            inner(v, f_u) * dx)

if include_viscosity:
    if include_viscosity_dgu:
        # Taken from http://maths.dur.ac.uk/~dma0mpj/summer_school/IPHO.pdf
        sigma = 1. # Penalty parameter.
        # Set taudgu=-1 for SIPG, taudgu=0 for IIPG, and taudgu=1 for NIPG
        taudgu = 1.
        edgelen = FacetArea(mesh)('+')  # Facetarea is continuous, so
        # we can select either side
        alpha = sigma/edgelen
        F_u_tent += nu * inner(grad(v), grad(u_mean)) * dx
        for d in range(2):
            F_u_tent += - nu * inner(avg(grad(u_mean[d])), jump(v[d], n))*dS
            F_u_tent += - nu * taudgu * inner(avg(grad(v[d])), jump(u[d], n))*dS
            F_u_tent += alpha * nu * inner(jump(u[d], n), jump(v[d], n))*dS
    else:
        F_u_tent += nu * inner(grad(v), grad(u_mean)) * dx


if include_advection:
    F_u_tent += inner(v, grad(u_mean)*u_bash) * dx

if include_supg:
    # tau1 = h/2/u_norm
    # tau1 = 0.5*h*pow(4.0*nu/h+2.0*u_norm, -1.0)
    tau1 = pow( (1./0.5/dt)**2 + (2.*u_norm/h)**2 + 9.0*(4.*nu/h/h)**2, -1.0/2)
    # tau1 = pow(1 / 0.5 / dt + 2. * u_norm / h + 4. * nu / h / h, -1.0)
    r1 = inner((1 / dt) * u_diff + g * grad(eta0) + friction / H * u_norm * u_mean - f_u, tau1 * grad(v) * u_bash) * dx
    if include_advection:
        r1 += inner(grad(u_mean) * u_bash, tau1 * grad(v) * u_bash) * dx
    if include_viscosity:
        r1 += inner(-div(nu * grad(u_mean)), tau1 * grad(v) * u_bash) * dx
    F_u_tent += r1


if include_ip:
    alpha = Constant(100.0)
    interior_penality = avg(alpha)*avg(h)**2*inner(jump(grad(u), n), jump(grad(v), n)) * dS
    F_u_tent += interior_penality


# Pressure correction step
eta_diff = eta - eta0
ut_mean = theta * ut + (1. - theta) * u0
F_p_corr = q*eta_diff*dx + g * dt**2 * theta**2 * H * inner(grad(q), grad(eta_diff))*dx - dt*H*inner(grad(q), ut_mean)*dx
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
u_file = File(prefix+"_u.pvd")
eta_file = File(prefix+"_eta.pvd")
eddy_vis = File(prefix+"_eddy.pvd")
u_file << (u1, float(t))
eta_file << (eta1, float(t))
eddy_vis << (project_nu, float(t))


while float(t - finish_time) <= - 1e3*DOLFIN_EPS:

    timestep += 1
    t = Constant(t+dt)
    print timestep

    # update source term
    t_theta = Constant(t - (1.0 - theta) * dt)

    # update boundary condition
    opentidal.t = t

    # update eddy_viscosity
    if include_les:
        log(INFO, "Compute eddy viscosity.")

        les.solve()
        eddy_viscosity = les.eddy_viscosity
        nu = Constant(viscosity) + eddy_viscosity
    project_nu = project(nu, Q)


    # Compute tentative velocity step
    log(INFO, "Solve for tentative velocity.")
    if include_crosswind:
        # tau2 = h/2/u_norm
        # tau2 = 0.5*h*pow(4.0*nu/h+2.0*u_norm, -1.0)
        tau2 = pow( (1./0.5/dt)**2 + (2.*u_norm/h)**2 + 9.0*(4.*nu/h/h)**2, -1.0/2)
        # tau2 = pow(1 / 0.5 / dt + 2. * u_norm / h + 4. * nu / h / h, -1.0)
        if assemble(du_norm*dx) > 1e-3:
            u_bash_parallel = inner(u_bash, grad(u0_mean[0])) * grad(u0_mean[0]) / du_norm ** 2
            r2u = inner((1 / dt) * u_diff[0] + g * grad(eta0)[0] + friction / H * u_norm * u_mean[0] - f_u[0],
                        tau2 * dot(grad(v[0]), u_bash_parallel)) * dx
            if include_advection:
                r2u += inner(dot(grad(u_mean[0]), u_bash), tau2 * dot(grad(v[0]), u_bash_parallel)) * dx
            if include_viscosity:
                r2u += inner(-div(nu * grad(u_mean[0])), tau2 * dot(grad(v[0]), u_bash_parallel)) * dx
            F_u_tent += r2u
        if assemble(dv_norm*dx) > 1e-3:
            v_bash_parallel = inner(u_bash, grad(u0_mean[1])) * grad(u0_mean[1]) / dv_norm ** 2
            r2v = inner((1 / dt) * u_diff[1] + g * grad(eta0)[1] + friction / H * u_norm * u_mean[1] - f_u[1],
                        tau2 * dot(grad(v[1]), v_bash_parallel)) * dx
            if include_advection:
                r2v += inner(dot(grad(u_mean[1]), u_bash), tau2 * dot(grad(v[1]), v_bash_parallel)) * dx
            if include_viscosity:
                r2v += inner(-div(nu * grad(u_mean[1])), tau2 * dot(grad(v[1]), v_bash_parallel)) * dx
            F_u_tent += r2v

    a_u_tent = lhs(F_u_tent)
    L_u_tent = rhs(F_u_tent)
    A_u_tent = assemble(a_u_tent)
    b = assemble(L_u_tent)
    for bc in bcu: bc.apply(A_u_tent, b)
    solve(A_u_tent, ut.vector(), b)

    if include_crosswind:
        if assemble(du_norm*dx) > 1e-3:
            F_u_tent -= r2u
        if assemble(dv_norm*dx) > 1e-3:
            F_u_tent -= r2v



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


    # Write to fileu_file << (u1, float(t))
    u_file << (u1, float(t))
    eta_file << (eta1, float(t))
    eddy_vis << (project_nu, float(t))


    # Rotate functions for next time_step
    u00.assign(u0)
    u0.assign(u1)
    eta0.assign(eta1)


log(INFO, "End of time loop.")
print "time_step is ", timestep
