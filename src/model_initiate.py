

from fenics import *
from netCDF4 import Dataset
from make_sto_basis import make_sto_basis
from make_sto_ijk import make_sto_ijk
from make_sto_modes import make_sto_modes
from make_mesh import make_mesh
from make_bath import make_bath
from make_ic import make_initial_object_list
from make_bc import make_boundary_object_list
from make_wind import MakeWind
from make_les import LES


parameters["form_compiler"]["cpp_optimize"] = True
ffc_options = {"optimize": True, "eliminate_zeros": True, "precompute_basis_const": True, "precompute_ip_const": True}


class ModelInitiate:
    """
    class to import mesh, set up initial condition and boundary condition, prepare wind input at each time step.

    Parameters:
    -----------
    inputs: Object
        object of class ModelReadInput from model_read_input.py.

    theta: float
        explicit implicit parameter.

    n_time_step : int
        number of time step.

    n_modes : int
        number of stochastic modes, it equals combination(sto_deg, sto_deg + sto_dim)

    B : FunctionSpace object
        object of Fenics.FunctionSpace class, used to build up shape function w.r.t. the domain.
        B is scalar functionspace in physical domain.

    C : FunctionSpace object
        object of Fenics.FunctionSpace class, used to build up shape function w.r.t. the domain.
        C is 2D vector functionspace in physical domain.

    W : FunctionSpace object
        object of Fenics.FunctionSpace class, used to build up shape function w.r.t. the domain.
        W is n_modes vector functionspace in stochastic domain.

    P : FunctionSpace object
        object of Fenics.FunctionSpace class, used to build up shape function w.r.t. the domain.
        B is 2 * n_modes vector functionspace in stochastic domain.

    v : Function object
        test function for water velocity of 2 * n_modes vector functionspace in stochastic domain

    u : Function object
        trial function of 2 * n_modes vector functionspace in stochastic domain

    q : Function object
        test function for surface elevation of n_modes vector functionspace in stochastic domain

    eta : Function object
        trial function of n_modes vector functionspace in stochastic domain.

    u00 : Function object
        water velocity at n-1 time step.

    u0 : Function object
        water velocity at n time step.

    ut : Function object
        water velocity at the tentative time step, which is between (n,n+1) time step.

    u1 : Function object
        water velocity at the (n+1)th time step.

    eta0 : Function object
        surface elevation at the nth time step.

    eta1 : Function object
        surface elevation at the (n+1)th time step.

    H : Function object
        total water column height.

    u_bc_object_list : list
        a list of velocity related boundary condition. like free slip, impenetration, etc.

    eta_bc_object_list : list
        a list of surface elevation boundary condition. like free slip, impenetration, etc.

    uTimeDependent : bool
        True only if any boundary conditions in the u_bc_object_list have time dependence.

    etaTimeDependent : bool
        True only if any boundary conditions in the eta_bc_object_list have time dependence.

    u_list_expression : object
        object of Fenics.Expression class, for the built of water velocity boundary condition.

    eta_list_expression : object
        object of FEnics.Expression class, for the built of surface elevation boundary condition.

    """

    def __init__(self, inputs):

        self.inputs = inputs

        self.theta = Constant(self.inputs.theta)
        self.g = Constant(self.inputs.Gravity)

        self.dt = Constant(self.inputs.time_step)
        self.finish_time = Constant(self.inputs.end_time)
        self.t = Constant(self.inputs.start_time)
        self.n_time_step = int((self.inputs.end_time - self.inputs.start_time) / self.inputs.time_step)

        self.basis_str = None
        self.n_modes = None
        self.ort_pol = None
        self.stoIJK = None

        self.mesh = None

        self.B = None
        self.C = None
        self.W = None
        self.P = None

        self.v = None
        self.u = None
        self.q = None
        self.eta = None

        self.u00 = None
        self.u0 = None
        self.ut = None
        self.u1 = None
        self.eta0 = None
        self.eta1 = None
        self.write_u1 = None
        self.write_eta1 = None

        self.H = None

        self.u_bc_object_list = None
        self.eta_bc_object_list = None
        self.uTimeDependent = None
        self.etaTimeDependent = None
        self.u_list_expression = None
        self.eta_list_expression = None

        self.wind = None
        self.wind_para_x = None
        self.wind_para_y = None
        self.pressure = None
        self.x_coord = None
        self.y_coord = None
        self.x_deg = None
        self.y_deg = None
        self.wind_vector = None

        self.windDrag_expression_object_list = []
        self.bottomDrag_expression_object_list = []
        self.nu_expression_object_list = []
        self.les = []

    def make_stochastic_basis(self):
        """ build up stochastic basis function in random space. """
        self.basis_str = make_sto_basis(self.inputs.dist_name, self.inputs.sto_poly_deg,
                                        self.inputs.sto_poly_dim, self.inputs.coefficient)
        self.stoIJK = make_sto_ijk(self.basis_str)
        # n_modes, ort_pol will be used by this class methods.
        self.n_modes = self.basis_str.get("n_modes")
        self.ort_pol = self.basis_str.get("basis")

    def make_mesh(self):
        """ build up mesh object which is compatible with Fenics."""
        self.mesh = make_mesh(self.inputs.domain)

    def make_function_space(self):
        """ prepare necessary function space for entire simulation. """
        vv = VectorElement(family='CG', cell=self.mesh.ufl_cell(), degree=2, dim=2)
        qq = FiniteElement(family='CG', cell=self.mesh.ufl_cell(), degree=1)
        self.B = FunctionSpace(self.mesh, qq)
        self.C = FunctionSpace(self.mesh, vv)
        string_v = "["
        string_q = "["
        for _ in range(self.n_modes - 1):
            string_v = string_v + "vv,"
            string_q = string_q + "qq,"
        string_v = string_v + "vv]"
        string_q = string_q + "qq]"
        if self.n_modes == 1:
            self.W = self.C
            self.P = self.B
        else:
            self.W = FunctionSpace(self.mesh, MixedElement(eval(string_v)))
            self.P = FunctionSpace(self.mesh, MixedElement(eval(string_q)))

    def make_function(self):
        """ Initialize input & solution functions for entire simulation. """
        self.v = TestFunction(self.W)
        self.u = TrialFunction(self.W)
        self.q = TestFunction(self.P)
        self.eta = TrialFunction(self.P)
        self.u00 = Function(self.W, name="u00")
        self.u0 = Function(self.W, name="u0")
        self.ut = Function(self.W, name="ut")
        self.u1 = Function(self.W, name="u")
        self.eta0 = Function(self.P, name="eta0")
        self.eta1 = Function(self.P, name="eta")
        # self.write_u1 = Function(self.C, name="u_modes")
        # self.write_eta1 = Function(self.B, name="eta_modes")
        self.wind_vector = Function(self.C, name="wind_xy")

    def make_bath(self):
        """ build up bathymetry object which is compatible with Fenics. """
        project_bottom_function = make_bath(self.inputs.bathymetry, self.basis_str, self.P)
        if self.inputs.linear_divergence:
            self.H = project_bottom_function
        else:
            self.H = project_bottom_function + self.eta0

    def make_ic(self):
        """ build up initial condition object which is compatible with Fenics. """
        if self.inputs.include_supg or self.inputs.include_crosswind:
            u_ic_function, eta_ic_function = make_initial_object_list(self.inputs.initial_u, self.inputs.initial_eta,
                                                                      self.basis_str, self.W, self.P, True)
        else:
            u_ic_function, eta_ic_function = make_initial_object_list(self.inputs.initial_u, self.inputs.initial_eta,
                                                                      self.basis_str, self.W, self.P, False)
        self.u0.assign(u_ic_function)
        self.u00.assign(u_ic_function)
        self.u1.assign(u_ic_function)
        self.eta0.assign(eta_ic_function)
        self.eta1.assign(eta_ic_function)

    def make_bc(self):
        """ build up boundary condition which is compatible with Fenics. """
        self.u_bc_object_list, self.eta_bc_object_list, self.uTimeDependent, self.etaTimeDependent, \
            self.u_list_expression, self.eta_list_expression = make_boundary_object_list(
                self.inputs.boundary_u, self.inputs.boundary_eta, self.basis_str, self.inputs.bc_file, self.mesh,
                self.inputs.domain, self.W, self.P, self.t, self.inputs.tidal_amplitude, self.inputs.tidal_period)

    def make_wind(self):
        """ build up 10 meter wind velocity object. """
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
        self.wind_para_x = Expression(cppcode=code, element=self.B.ufl_element(), domain=self.mesh)
        self.wind_para_y = Expression(cppcode=code, element=self.B.ufl_element(), domain=self.mesh)
        self.pressure = Expression(cppcode=code, element=self.B.ufl_element(), domain=self.mesh)
        nc = Dataset(self.inputs.input_dir + self.inputs.bath_file, 'r', format='NETCDF4')
        self.x_coord = nc.variables['x_coord'][:].data.tolist()
        self.y_coord = nc.variables['y_coord'][:].data.tolist()
        self.x_deg = nc.variables['lon'][:].data.tolist()
        self.y_deg = nc.variables['lat'][:].data.tolist()
        for st in range(len(self.x_coord)):
            self.wind_para_x.initial(self.x_coord[st], self.y_coord[st], 0.0)
            self.wind_para_y.initial(self.x_coord[st], self.y_coord[st], 0.0)
            self.pressure.initial(self.x_coord[st], self.y_coord[st], 1013.0)
        # self.wind = MakeWind(self.inputs.input_dir + self.inputs.wind_file, self.inputs.wind_dt)
        self.wind = MakeWind(self.inputs)

    def make_const_wind(self):
        self.wind_para_x = Constant(self.inputs.wind_x)
        self.wind_para_y = Constant(self.inputs.wind_y)

    def make_stochastic_parameter(self):
        """ build up entire stochastic modes for uncertain inputs. """
        nu_expression_list, bottom_drag_expression_list, wind_drag_expression_list = make_sto_modes(
            self.basis_str, self.inputs.sto_viscosity, self.inputs.sto_bottomDrag, self.inputs.sto_windDrag)

        if self.inputs.include_wind_stress:
            for kk in range(self.n_modes):
                self.windDrag_expression_object_list.append(
                    Expression(wind_drag_expression_list[kk], element=self.B.ufl_element(), domain=self.mesh))

        for kk in range(self.n_modes):
            self.bottomDrag_expression_object_list.append(
                Expression(bottom_drag_expression_list[kk], element=self.B.ufl_element(), domain=self.mesh))

        if self.inputs.include_les:
            if self.n_modes == 1:
                self.les.append(LES(self.B, self.u0, self.inputs.les_parameters['smagorinsky_coefficient']))
            else:
                for kk in range(self.n_modes):
                    self.les.append(LES(self.B, self.u0.split()[kk],
                                        self.inputs.les_parameters['smagorinsky_coefficient']))
            for kk in range(self.n_modes):
                self.nu_expression_object_list.append(Expression(nu_expression_list[kk], element=self.B.ufl_element(),
                                                                 domain=self.mesh) + self.les[kk].eddy_viscosity)
        else:
            for kk in range(self.n_modes):
                self.nu_expression_object_list.append(Expression(nu_expression_list[kk], element=self.B.ufl_element(),
                                                                 domain=self.mesh))
