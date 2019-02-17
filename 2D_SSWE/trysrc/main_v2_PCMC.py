
import sys
sys.path.insert(0, '/workspace/Documentation/Research_Doc/SFEM_Doc/7-NS-github/input')
from input_generalized import *
from makestobasis import MakeStoBasis
from makestomodes import MakeStoModes
from makestoijk import MakeStoIJK
from makemesh import MakeMesh
from make_bath import make_bath
from makeic import MakeInitialObjectFunc
from makebc import MakeBoundaryObjectList
from makewind import MakeWind
from netCDF4 import Dataset
from fenics import *
import numpy as np


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
for i in range(nmodes - 1):
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
nu_expression_list, bottomDrag_expression_list, windDrag_expression_list = MakeStoModes(basis_struct, sto_viscosity,
                                                                                        sto_bottomDrag, sto_windDrag)

# Get time parameter setting
theta = Constant(theta)
dt = Constant(timestep)
finish_time = Constant(endtime)
t = Constant(startime)
n_time_step = int((endtime - startime) / timestep)

# Get equation parameter setting
g = Constant(Gravity)

# Get wind parameter setting
B = FunctionSpace(mesh, Q)
C = VectorFunctionSpace(mesh, "CG", 2, dim=2)
if include_wind:
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

            void initial(double _x, double _y, double _v)
            {
                node *newNode = new node;
                newNode->x = _x;
                newNode->y = _y;
                newNode->v = _v;
                newNode->next = root;
                root = newNode;
                //cout << _v << endl;
            };

            void update(double _x, double _y, double _v)
            {
                node *p = findp(_x, _y);
                //cout << p->v << "  " << _v << endl;
                p->v = _v;
            };

            void eval(Array<double>& values, const Array<double>& x) const
            {
                double vv = findval(x[0], x[1]);
                values[0] = vv;
                //cout << x[0] << "  " << x[1] << "  " << vv << endl;
                //cout << "  " << endl;
            };

            node * findp(double _x, double _y) const
            {
                node *p = root;
                while (p->next != 0)
                {
                    if ( (fabs(p->x - _x)<1.0e-4) && (fabs(p->y - _y)<1.0e-4) )
                    {
                        return p;
                    }
                    else
                    {
                        p = p->next;
                    }
                }
                return 0;
            }

            double findval(double _x, double _y) const
            {
                node *p = root;
                while (p->next != 0)
                {   
                    if ( (fabs(p->x - _x)<1.0e-4) && (fabs(p->y - _y)<1.0e-4) )
                    {
                        //cout << fabs(p->x-_x) << "  " << fabs(p->y-_y) << endl;
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
    wind_para_x = Expression(cppcode=code, element=B.ufl_element())
    wind_para_y = Expression(cppcode=code, element=B.ufl_element())
    nc = Dataset(input_dir + bathyfile, 'r', format='NETCDF4')
    x_coord = nc.variables['x_coord'][:].data.tolist()
    y_coord = nc.variables['y_coord'][:].data.tolist()
    x_deg = nc.variables['lon'][:].data.tolist()
    y_deg = nc.variables['lat'][:].data.tolist()
    for st in range(len(x_coord)):
        wind_para_x.initial(x_coord[st], y_coord[st], 0.0)
        wind_para_y.initial(x_coord[st], y_coord[st], 0.0)
    wind = MakeWind(input_dir + windfile, time_step=wind_dt)

# Get the bathymetry
project_bottom_function = make_bath(bathymetry, basis_struct, P)

# plot(project_bottom_function)

# Get initial condition in stochastic setting
u_ic_function, eta_ic_function = MakeInitialObjectFunc(initial_u, initial_eta, basis_struct, W, P)

# Get boundary condition in stochastic setting
u_bc_object_list, eta_bc_object_list, uTimeDependent, etaTimeDependent, u_list_expression, eta_list_expression = MakeBoundaryObjectList(
    boundary_u, boundary_eta, basis_struct, bcfile, mesh, domain, W, P, t, project_bottom_function)


##################################################################
# starting initializing
##################################################################

# For large eddy model ==> error: grad(u) + grad(u).
# T shape is different, can't match at all. ==> consider: build a new subspace.
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
ut = Function(W)  # Tentative velocity
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
windDrag_expression_object_list = []
bottomDrag_expression_object_list = []
nu_expression_object_list = []
les = []
for k in range(nmodes):
    if include_wind:
        windDrag_expression_object_list.append(Expression(windDrag_expression_list[k], element=B.ufl_element()))

    bottomDrag_expression_object_list.append(Expression(bottomDrag_expression_list[k], element=B.ufl_element()))

    if include_les:
        if nmodes == 1:
            les.append(LES(B, u0, les_parameters['smagorinsky_coefficient']))
        else:
            les.append(LES(B, u0.split()[k], les_parameters['smagorinsky_coefficient']))
        nu_expression_object_list.append(
            Expression(nu_expression_list[k], element=B.ufl_element()) + les[k].eddy_viscosity)
    else:
        nu_expression_object_list.append(Expression(nu_expression_list[k], element=B.ufl_element()))

# Tentative Velocity step
u_mean = theta * u + (1. - theta) * u0
u_bash = 3. / 2 * u0 - 1. / 2 * u00
u_diff = u - u0
norm_u0 = inner(u0, u0) ** 0.5
F_u_tent = 0.0

for k in range(nmodes):
    if nmodes == 1:
        F_u_tent = F_u_tent + (1 / dt) * inner(v, u_diff) * dx + \
                   g * inner(v, grad(eta0)) * dx
    else:
        F_u_tent = F_u_tent + (1 / dt) * inner(v[2 * k], u_diff[2 * k]) * dx + \
                   (1 / dt) * inner(v[2 * k + 1], u_diff[2 * k + 1]) * dx + \
                   g * inner(v[2 * k], grad(eta0[k])[0]) * dx + \
                   g * inner(v[2 * k + 1], grad(eta0[k])[1]) * dx
    for j in range(nmodes):
        for i in range(nmodes):
            if abs(stoIJK[i][j][k]) > DOLFIN_EPS:
                if nmodes == 1:
                    F_u_tent = F_u_tent + (
                            stoIJK[i][j][k] * inner(u0, u0) ** 0.5 / H * bottomDrag_expression_object_list[
                        i] * inner(u_mean[2 * j], v[2 * k]) * dx +
                            stoIJK[i][j][k] * inner(u0, u0) ** 0.5 / H * bottomDrag_expression_object_list[
                                i] * inner(u_mean[2 * j + 1], v[2 * k + 1]) * dx)
                else:
                    F_u_tent = F_u_tent + (stoIJK[i][j][k] * (u0[0] ** 2 + u0[1] ** 2) ** 0.5 / H[0] *
                                           bottomDrag_expression_object_list[i] * inner(u_mean[2 * j], v[2 * k]) * dx +
                                           stoIJK[i][j][k] * (u0[0] ** 2 + u0[1] ** 2) ** 0.5 / H[0] *
                                           bottomDrag_expression_object_list[i] * inner(u_mean[2 * j + 1],
                                                                                        v[2 * k + 1]) * dx)
                if include_viscosity:
                    F_u_tent = F_u_tent + (stoIJK[i][j][k] * nu_expression_object_list[i] * inner(grad(v[2 * k]), grad(
                        u_mean[2 * j])) * dx +
                                           stoIJK[i][j][k] * nu_expression_object_list[i] * inner(grad(v[2 * k + 1]),
                                                                                                  grad(u_mean[
                                                                                                           2 * j + 1])) * dx)
                if include_advection:
                    F_u_tent = F_u_tent + (stoIJK[i][j][k] * inner(v[2 * k],
                                                                   u_bash[2 * i] * grad(u_mean[2 * j])[0] + u_bash[
                                                                       2 * i + 1] * grad(u_mean[2 * j])[1]) * dx +
                                           stoIJK[i][j][k] * inner(v[2 * k + 1],
                                                                   u_bash[2 * i] * grad(u_mean[2 * j + 1])[0] + u_bash[
                                                                       2 * i + 1] * grad(u_mean[2 * j + 1])[1]) * dx)
    if include_wind:

        norm_wind = (wind_para_x ** 2 + wind_para_y ** 2) ** 0.5

        if nmodes == 1:
            F_u_tent = F_u_tent - v[2 * k] * rho_air / rho_water * norm_wind * wind_para_x / H * windDrag_expression_object_list[k] * dx
            F_u_tent = F_u_tent - v[2 * k + 1] * rho_air / rho_water * norm_wind * wind_para_y / H * windDrag_expression_object_list[k] * dx
        else:
            F_u_tent = F_u_tent - v[2 * k] * rho_air / rho_water * norm_wind * wind_para_x / H[0] * windDrag_expression_object_list[k] * dx
            F_u_tent = F_u_tent - v[2 * k + 1] * rho_air / rho_water * norm_wind * wind_para_y / H[0] * windDrag_expression_object_list[k] * dx

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
            if abs(stoIJK[i][j][k]) > DOLFIN_EPS:
                if nmodes == 1:
                    F_p_corr = F_p_corr + stoIJK[i][j][k] * theta ** 2 * g * dt ** 2 * H * inner(grad(q),
                                                                                                 grad(eta_diff)) * dx \
                               - stoIJK[i][j][k] * dt * H * inner(grad(q), ut_mean) * dx  # \
                    # + 1./h**sigma*inner(dot(ut_mean, n), dot())
                    # + stoIJK[i][j][k] * dt * q * grad(H*ut_mean[0])[0] * dx
                    # + stoIJK[i][j][k] * dt * q * grad(H*ut_mean[1])[1] * dx)
                else:
                    F_p_corr = F_p_corr + stoIJK[i][j][k] * theta ** 2 * g * dt ** 2 * H[i] * inner(grad(q[k]), grad(
                        eta_diff[j])) * dx \
                               - stoIJK[i][j][k] * dt * H[i] * grad(q[k])[0] * ut_mean[2 * j] * dx \
                               - stoIJK[i][j][k] * dt * H[i] * grad(q[k])[1] * ut_mean[2 * j + 1] * dx
                    # + stoIJK[i][j][k] * dt * q[k] * grad(H[i]*ut_mean[2*j])[0] * dx \
                    # + stoIJK[i][j][k] * dt * q[k] * grad(H[i]*ut_mean[2*j+1])[1] * dx
a_p_corr = lhs(F_p_corr)
L_p_corr = rhs(F_p_corr)

# Velocity correction step
eta_diff = eta1 - eta0
F_u_corr = 0.0
if nmodes == 1:
    F_u_corr = inner(v, u) * dx - inner(v, ut) * dx + dt * g * theta * inner(v, grad(eta_diff)) * dx
else:
    for k in range(nmodes):
        F_u_corr = F_u_corr + (
                u[2 * k] * v[2 * k] - ut[2 * k] * v[2 * k] + dt * g * theta * v[2 * k] * grad(eta_diff[k])[0]) * dx \
                   + (u[2 * k + 1] * v[2 * k + 1] - ut[2 * k + 1] * v[2 * k + 1] + dt * g * theta * v[2 * k + 1] *
                      grad(eta_diff[k])[1]) * dx
a_u_corr = lhs(F_u_corr)
L_u_corr = rhs(F_u_corr)

# Assemble matrices
a_u_tent_solver = KrylovSolver("gmres", "ilu")
a_u_tent_solver.parameters['absolute_tolerance'] = 1E-4
a_u_tent_solver.parameters['relative_tolerance'] = 1E-3
a_u_tent_solver.parameters['maximum_iterations'] = 1000
a_u_tent_solver.parameters['monitor_convergence'] = True
a_u_tent_solver.parameters['nonzero_initial_guess'] = True

if linear_divergence:
    A_p_corr = assemble(a_p_corr)
    for bc in eta_bc_object_list: bc.apply(A_p_corr)
    a_p_corr_solver = LUSolver(A_p_corr)
    a_p_corr_solver.parameters["reuse_factorization"] = True
else:
    a_p_corr_solver = KrylovSolver("gmres", "ilu")
    a_p_corr_solver.parameters['absolute_tolerance'] = 1E-4
    a_p_corr_solver.parameters['relative_tolerance'] = 1E-3
    a_p_corr_solver.parameters['maximum_iterations'] = 1000
    a_p_corr_solver.parameters['monitor_convergence'] = True

A_u_corr = assemble(a_u_corr)
for bc in u_bc_object_list: bc.apply(A_u_corr)
# a_u_corr_solver = KrylovSolver(A_u_corr, "gmres", "ilu")
# a_u_corr_solver.parameters['absolute_tolerance'] = 1E-7
# a_u_corr_solver.parameters['relative_tolerance'] = 1E-4
# a_u_corr_solver.parameters['maximum_iterations'] = 1000
# a_u_corr_solver.parameters['monitor_convergence'] = True
a_u_corr_solver = LUSolver(A_u_corr)
a_u_corr_solver.parameters["reuse_factorization"] = True

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
w = Function(C, name="wind_speed")

if nmodes == 1:
    if USEpvd:
        u_file = File(outputdir + "u_" + outputstr + "only_mode_0.pvd")
        eta_file = File(outputdir + "eta_" + outputstr + "only_mode_0.pvd")
        # eddy_file = File(output_dir + "eddy_" + outputstr + "only_mode_0.pvd")
        u_file << (u1, float(t))
        eta_file << (eta1, float(t))
    else:
        u_file = XDMFFile(outputdir + "u_" + outputstr + "only_mode_0.xdmf")
        eta_file = XDMFFile(outputdir + "eta_" + outputstr + "only_mode_0.xdmf")
        u_file.write(u1, float(t))
        eta_file.write(eta1, float(t))

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
            u_file.append(XDMFFile(outputdir + "u_" + outputstr + "{:02d}".format(mode) + ".xdmf"))
            eta_file.append(XDMFFile(outputdir + "eta_" + outputstr + "{:02d}".format(mode) + ".xdmf"))

        tempu1.assign(u1splits[mode])
        tempeta1.assign(eta1splits[mode])

        if USEpvd:
            u_file[mode] << (tempu1, float(t))
            eta_file[mode] << (tempeta1, float(t))
        else:
            u_file[mode].write(tempu1, float(t))
            eta_file[mode].write(tempeta1, float(t))

if include_wind:
    if USEpvd:
        wind_xy_output_file = File(outputdir + "wvel_" + outputstr + "0.pvd")
    else:
        wind_xy_output_file = XDMFFile(outputdir + "wvel_" + outputstr + "0.xdmf")


# Prepare binrandom for PC
testnodes = [[a, b] for a in testnodeX for b in testnodeY]
samplex = np.linspace(coefficient[0], coefficient[1], nsample)
sampley = np.linspace(coefficient[2], coefficient[3], nsample)
sampleX, sampleY = np.meshgrid(samplex, sampley)
binrandom_u1 = np.zeros([nsample, nsample, n_time_step + 1, len(testnodes)])
binrandom_v1 = np.zeros([nsample, nsample, n_time_step + 1, len(testnodes)])
binrandom_eta1 = np.zeros([nsample, nsample, n_time_step + 1, len(testnodes)])
orthpol = basis_struct.get("basis")

for j in range(nsample):
    for k in range(nsample):
        orth_list = [orthpol[mode](samplex[j], sampley[k]) for mode in range(nmodes)]
        for m in range(len(testnodes)):
            u1_list = [u1(testnodes[m][0], testnodes[m][1])[2 * p] for p in range(nmodes)]
            v1_list = [u1(testnodes[m][0], testnodes[m][1])[2 * p + 1] for p in range(nmodes)]
            binrandom_u1[j, k, timestepcount, m] = np.dot(orth_list, u1_list)
            binrandom_v1[j, k, timestepcount, m] = np.dot(orth_list, v1_list)
            binrandom_eta1[j, k, timestepcount, m] = np.dot(orth_list, eta1(testnodes[m][0], testnodes[m][1]))

# Start timestepping.
n = FacetNormal(mesh)
total_mass = []
time_stamp = []

while float(t - finish_time) < -1e-3:

    time_stamp.append(float(t))

    # update current time t.==> successfully updated
    timestepcount += 1
    t = Constant(t + dt)
    print timestepcount
    print "time is: ", float(t)

    # update wind force term with time t.
    t_theta = Constant(t - (1 - theta) * dt)
    if include_wind:
        wind.get_prepared(current_time=float(t))
        wind_para_x_list = []
        wind_para_y_list = []
        for ts in range(len(x_deg)):
            wind_x, wind_y, pressure = wind.get_wind_x_wind_y(x_coord_deg=x_deg[ts], y_coord_deg=y_deg[ts])
            wind_para_x_list.append(wind_x)
            wind_para_y_list.append(wind_y)
        for sm in range(len(x_coord)):
            wind_para_x.update(x_coord[sm], y_coord[sm], wind_para_x_list[sm])
            wind_para_y.update(x_coord[sm], y_coord[sm], wind_para_y_list[sm])
        temp = project(as_vector([wind_para_x, wind_para_y]), C)
        w.assign(temp)
        wind_xy_output_file << (w, float(t))

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
            # nu_expression_object_list[k] = Expression(nu_expression_list[k], element=B.ufl_element()) + les[k].eddy_viscosity

            if DEBUG_mode:
                cc = project(nu_expression_object_list[k], B)
                plot(cc)

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

    # a_u_tent_solver.set_operator(A_u_tent)
    # a_u_tent_solver.solve(ut.vector(), b)
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

        # a_p_corr_solver.set_operator(A_p_corr)
        # a_p_corr_solver.solve(eta1.vector(), b)
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

    print("===========================================")

    # Write to file
    if nmodes == 1:
        if USEpvd:
            u_file << (u1, float(t))
            eta_file << (eta1, float(t))
        else:
            u_file.write(u1, float(t))
            eta_file.write(eta1, float(t))
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
    for j in range(nsample):
        for k in range(nsample):
            orth_list = [orthpol[mode](samplex[j], sampley[k]) for mode in range(nmodes)]
            for m in range(len(testnodes)):
                u1_list = [u1(testnodes[m][0], testnodes[m][1])[2 * p] for p in range(nmodes)]
                v1_list = [u1(testnodes[m][0], testnodes[m][1])[2 * p + 1] for p in range(nmodes)]
                binrandom_u1[j, k, timestepcount, m] = np.dot(orth_list, u1_list)
                binrandom_v1[j, k, timestepcount, m] = np.dot(orth_list, v1_list)
                binrandom_eta1[j, k, timestepcount, m] = np.dot(orth_list, eta1(testnodes[m][0], testnodes[m][1]))
            print "PC finish: " + str(j) + str(k)

    # check mass conservation.
    total_mass.append(assemble(H*dx))

time_stamp.append(float(t))
total_mass.append(assemble(H*dx))

log(INFO, "End of time loop.")
print "total time_step is ", timestepcount
# print "total mass is: ". total_mass

np.save(outputdir + "binrandom_eta1_all_points_order_" + str(sto_poly_deg), binrandom_eta1)
np.save(outputdir + "binrandom_u1_all_points_order_" + str(sto_poly_deg), binrandom_u1)
np.save(outputdir + "binrandom_v1_all_points_order_" + str(sto_poly_deg), binrandom_v1)
np.save(outputdir + "total_mass_at_every_time_step", total_mass)
np.save(outputdir + "time_stamp_at_every_time_step", time_stamp)
