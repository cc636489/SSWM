

from fenics import inner, grad, DOLFIN_EPS, dx, CellSize, div, dot


class DetModelFormulate:
    """
    Deterministic shallow water equation weak form formulation.

    Parameters:
    -----------
    inputs : object
        object of class ModelReadInput from model_read_input.py.

    initiate : object
        object of class ModelInitiate from model_initiate.py.

    u_mean : function object

    u_bash : function object

    u_diff : function object

    u0_norm : function object
        vector norm of u0 at the nth time step.

    eta_diff : function object

    ut_mean : function object

    eta1_diff : function object

    norm_wind : function object
        vector norm of wind velocity

    F_u_tent : form object
        weak form for solving tentative velocity

    F_p_corr : form object
        weak form for solving pressure

    F_u_corr : form object
        weak form for solving the (n+1)th water velocity

    h : object
        class object to give each element size

    tau1 : object
        standard SUPG stabilized parameter

    tau2 : object
        cross wind diffusivity stabilized parameter

    """

    def __init__(self, inputs, initiate):

        self.inputs = inputs
        self.initiate = initiate

        self.u_mean = self.initiate.theta * self.initiate.u + (1. - self.initiate.theta) * self.initiate.u0
        self.u_bash = 3. / 2 * self.initiate.u0 - 1. / 2 * self.initiate.u00
        self.u_diff = self.initiate.u - self.initiate.u0
        self.u0_norm = inner(self.initiate.u0, self.initiate.u0) ** 0.5

        self.eta_diff = self.initiate.eta - self.initiate.eta0
        self.ut_mean = self.initiate.theta * self.initiate.ut + (1. - self.initiate.theta) * self.initiate.u0

        self.eta1_diff = self.initiate.eta1 - self.initiate.eta0

        self.norm_wind = (self.initiate.wind_para_x ** 2 + self.initiate.wind_para_y ** 2) ** 0.5

        self.F_u_tent = 0
        self.F_p_corr = 0
        self.F_u_corr = 0

        self.h = CellSize(self.initiate.mesh)
        # self.tau1 = 0.5*h*pow(4.0*nu/h+2.0*u_norm, -1.0)
        # self.tau1 = pow( (1./0.5/dt)**2 + (2.*u_norm/h)**2 + 9.0*(4.*nu/h/h)**2, -1.0/2)
        self.tau1 = pow(1 / 0.5 / self.initiate.dt + 2. * self.u0_norm / self.h +
                        4. * self.initiate.nu_expression_object_list[0] / self.h / self.h, -1.0)
        self.tau2 = self.tau1

        self.F_u_tent += ((1 / self.initiate.dt) * inner(self.initiate.v, self.u_diff) * dx +
                          self.initiate.g * inner(self.initiate.v, grad(self.initiate.eta0)) * dx)

        self.F_p_corr += (inner(self.eta_diff, self.initiate.q) * dx + self.initiate.theta ** 2 *
                          self.initiate.g * self.initiate.dt ** 2 * self.initiate.H *
                          inner(grad(self.initiate.q), grad(self.eta_diff)) * dx -
                          self.initiate.dt * self.initiate.H * inner(grad(self.initiate.q), self.ut_mean) * dx)

        self.F_u_corr += (inner(self.initiate.v, self.initiate.u) * dx -
                          inner(self.initiate.v, self.initiate.ut) * dx +
                          self.initiate.dt * self.initiate.g * self.initiate.theta *
                          inner(self.initiate.v, grad(self.eta1_diff)) * dx)

    def add_convection(self):
        """ add convective term to tentative water velocity F_u_tent weak form. """
        self.F_u_tent += inner(self.initiate.v, grad(self.u_mean) * self.u_bash) * dx

    def add_viscosity(self):
        """ add viscosity term to tentative water velocity F_u_tent weak form. """
        self.F_u_tent += self.initiate.nu_expression_object_list[0] * \
                         inner(grad(self.initiate.v), grad(self.u_mean)) * dx

    def add_bottom_stress(self):
        """ add bottom friction source term to tentative water velocity F_u_tent weak form. """
        self.F_u_tent += self.u0_norm / self.initiate.H * self.initiate.bottomDrag_expression_object_list[0] * \
                         inner(self.u_mean, self.initiate.v) * dx

    def add_wind_stress(self):
        """ add wind stress source term to tentative water velocity F_u_tent weak form. """
        self.F_u_tent -= (self.initiate.v[0] * self.inputs.rho_air / self.inputs.rho_water *
                          self.norm_wind * self.initiate.wind_para_x / self.initiate.H *
                          self.initiate.windDrag_expression_object_list[0] * (0.75 + 0.067 * self.norm_wind)) * dx
        self.F_u_tent -= (self.initiate.v[1] * self.inputs.rho_air / self.inputs.rho_water *
                          self.norm_wind * self.initiate.wind_para_y / self.initiate.H *
                          self.initiate.windDrag_expression_object_list[0] * (0.75 + 0.067 * self.norm_wind)) * dx

    def add_atmospheric_pressure(self):
        """ add atmospheric pressure source term to tentative water velocity F_u_tent weak form. """
        self.F_u_tent += inner(self.initiate.v, grad(self.initiate.pressure / self.inputs.rho_water)) * dx

    def add_su_pg(self):
        """ add standard stabilized term to tentative water velocity F_u_tent weak form. """
        r = inner((1 / self.initiate.dt) * self.u_diff + self.initiate.g * grad(self.initiate.eta0),
                  self.tau1 * grad(self.initiate.v) * self.u_bash) * dx

        if self.inputs.include_convection:
            r += inner(grad(self.u_mean) * self.u_bash, self.tau1 * grad(self.initiate.v) * self.u_bash) * dx

        if self.inputs.include_viscosity:
            r += inner(-div(self.initiate.nu_expression_object_list[0] * grad(self.u_mean)),
                       self.tau1 * grad(self.initiate.v) * self.u_bash) * dx

        if self.inputs.include_bottom_stress:
            r += inner(self.u0_norm / self.initiate.H * self.initiate.bottomDrag_expression_object_list[0] *
                       self.u_mean, self.tau1 * grad(self.initiate.v) * self.u_bash) * dx

        if self.inputs.include_wind_stress:
            r -= (self.tau1 * dot(self.u_bash, grad(self.initiate.v[0])) * self.inputs.rho_air /
                  self.inputs.rho_water * self.u0_norm * self.initiate.wind_para_x / self.initiate.H *
                  self.initiate.windDrag_expression_object_list[0] * (0.75 + 0.067 * self.norm_wind)) * dx
            r -= (self.tau1 * dot(self.u_bash, grad(self.initiate.v[1])) * self.inputs.rho_air /
                  self.inputs.rho_water * self.u0_norm * self.initiate.wind_para_y / self.initiate.H *
                  self.initiate.windDrag_expression_object_list[0] * (0.75 + 0.067 * self.norm_wind)) * dx

        if self.inputs.include_atmospheric_pressure:
            r += inner(grad(self.initiate.pressure / self.inputs.rho_water),
                       self.tau1 * grad(self.initiate.v) * self.u_bash) * dx

        self.F_u_tent += r

    def add_crosswind(self):
        """ add cross wind diffusion term to tentative water velocity F_u_tent weak form. """
        u0_mean = self.initiate.theta * self.initiate.u0 + (1. - self.initiate.theta) * self.initiate.u00
        du_norm = inner(grad(self.initiate.u0[0]), grad(self.initiate.u0[0])) ** 0.5
        dv_norm = inner(grad(self.initiate.u0[1]), grad(self.initiate.u0[1])) ** 0.5

        u_bash_parallel = inner(self.u_bash, grad(u0_mean[0])) * grad(u0_mean[0]) / du_norm ** 2
        v_bash_parallel = inner(self.u_bash, grad(u0_mean[1])) * grad(u0_mean[1]) / dv_norm ** 2

        ru = inner((1 / self.initiate.dt) * self.u_diff[0] + self.initiate.g * grad(self.initiate.eta0)[0],
                   self.tau2 * dot(grad(self.initiate.v[0]), u_bash_parallel)) * dx
        rv = inner((1 / self.initiate.dt) * self.u_diff[1] + self.initiate.g * grad(self.initiate.eta0)[1],
                   self.tau2 * dot(grad(self.initiate.v[1]), v_bash_parallel)) * dx

        if self.inputs.include_convection:
            ru += inner(dot(grad(self.u_mean[0]), self.u_bash), self.tau2 *
                        dot(grad(self.initiate.v[0]), u_bash_parallel)) * dx
            rv += inner(dot(grad(self.u_mean[1]), self.u_bash), self.tau2 *
                        dot(grad(self.initiate.v[1]), v_bash_parallel)) * dx

        if self.inputs.include_viscosity:
            ru += inner(-div(self.initiate.nu_expression_object_list[0] * grad(self.u_mean[0])),
                        self.tau2 * dot(grad(self.initiate.v[0]), u_bash_parallel)) * dx
            rv += inner(-div(self.initiate.nu_expression_object_list[0] * grad(self.u_mean[1])),
                        self.tau2 * dot(grad(self.initiate.v[1]), v_bash_parallel)) * dx

        if self.inputs.include_bottom_stress:
            ru += inner(self.u0_norm / self.initiate.H * self.initiate.bottomDrag_expression_object_list[0] *
                        self.u_mean[0], self.tau2 * dot(grad(self.initiate.v[0]), u_bash_parallel)) * dx
            rv += inner(self.u0_norm / self.initiate.H * self.initiate.bottomDrag_expression_object_list[0] *
                        self.u_mean[1], self.tau2 * dot(grad(self.initiate.v[1]), v_bash_parallel)) * dx

        if self.inputs.include_wind_stress:
            ru -= (self.tau2 * dot(grad(self.initiate.v[0]), u_bash_parallel) * self.inputs.rho_air /
                   self.inputs.rho_water * self.u0_norm * self.initiate.wind_para_x / self.initiate.H *
                   self.initiate.windDrag_expression_object_list[0] * (0.75 + 0.067 * self.norm_wind)) * dx
            rv -= (self.tau2 * dot(grad(self.initiate.v[1]), v_bash_parallel) * self.inputs.rho_air /
                   self.inputs.rho_water * self.u0_norm * self.initiate.wind_para_y / self.initiate.H *
                   self.initiate.windDrag_expression_object_list[0] * (0.75 + 0.067 * self.norm_wind)) * dx

        if self.inputs.include_atmospheric_pressure:
            ru += inner(grad(self.initiate.pressure / self.inputs.rho_water)[0],
                        self.tau2 * dot(grad(self.initiate.v[0]), u_bash_parallel)) * dx
            rv += inner(grad(self.initiate.pressure / self.inputs.rho_water)[1],
                        self.tau2 * dot(grad(self.initiate.v[1]), v_bash_parallel)) * dx

        self.F_u_tent += ru + rv


class StoModelFormulate:
    """
    Stochastic shallow water equation weak form formulation.

    Parameters:
    -----------
    inputs : object
        object of class ModelReadInput from model_read_input.py.

    initiate : object
        object of class ModelInitiate from model_initiate.py.

    u_mean : function object

    u_bash : function object

    u_diff : function object

    u0_norm : function object
        vector norm of u0 at the nth time step.

    eta_diff : function object

    ut_mean : function object

    eta1_diff : function object

    norm_wind : function object
        vector norm of wind velocity

    F_u_tent : form object
        weak form for solving tentative velocity

    F_p_corr : form object
        weak form for solving pressure

    F_u_corr : form object
        weak form for solving the (n+1)th water velocity

    h : object
        class object to give each element size

    tau1 : object
        standard SUPG stabilized parameter

    tau2 : object
        cross wind diffusivity stabilized parameter

    """

    def __init__(self, inputs, initiate):

        self.inputs = inputs
        self.initiate = initiate

        self.u_mean = self.initiate.theta * self.initiate.u + (1. - self.initiate.theta) * self.initiate.u0
        self.u_bash = 3. / 2 * self.initiate.u0 - 1. / 2 * self.initiate.u00
        self.u_diff = self.initiate.u - self.initiate.u0
        self.norm_u0 = inner(self.initiate.u0, self.initiate.u0) ** 0.5

        self.eta_diff = self.initiate.eta - self.initiate.eta0
        self.ut_mean = self.initiate.theta * self.initiate.ut + (1. - self.initiate.theta) * self.initiate.u0

        self.eta1_diff = self.initiate.eta1 - self.initiate.eta0

        self.norm_wind = (self.initiate.wind_para_x ** 2 + self.initiate.wind_para_y ** 2) ** 0.5

        self.F_u_tent = 0
        self.F_p_corr = 0
        self.F_u_corr = 0
        self.R = 0

        for k in range(self.initiate.n_modes):
            self.F_u_tent += (1 / self.initiate.dt) * inner(self.initiate.v[2 * k], self.u_diff[2 * k]) * dx + \
                             (1 / self.initiate.dt) * inner(self.initiate.v[2 * k + 1], self.u_diff[2 * k + 1]) * dx + \
                             self.initiate.g * inner(self.initiate.v[2 * k], grad(self.initiate.eta0[k])[0]) * dx + \
                             self.initiate.g * inner(self.initiate.v[2 * k + 1], grad(self.initiate.eta0[k])[1]) * dx

            self.F_p_corr += inner(self.eta_diff[k], self.initiate.q[k]) * dx

            self.F_u_corr += (self.initiate.u[2 * k] * self.initiate.v[2 * k] -
                              self.initiate.ut[2 * k] * self.initiate.v[2 * k] +
                              self.initiate.dt * self.initiate.g * self.initiate.theta *
                              self.initiate.v[2 * k] * grad(self.eta1_diff[k])[0]) * dx + \
                             (self.initiate.u[2 * k + 1] * self.initiate.v[2 * k + 1] -
                              self.initiate.ut[2 * k + 1] * self.initiate.v[2 * k + 1] +
                              self.initiate.dt * self.initiate.g * self.initiate.theta *
                              self.initiate.v[2 * k + 1] * grad(self.eta1_diff[k])[1]) * dx
            for j in range(self.initiate.n_modes):
                for i in range(self.initiate.n_modes):
                    if abs(self.initiate.stoIJK[i][j][k]) > DOLFIN_EPS:
                        self.F_p_corr += self.initiate.stoIJK[i][j][k] * self.initiate.theta ** 2 * \
                                         self.initiate.g * self.initiate.dt ** 2 * self.initiate.H[i] * \
                                         inner(grad(self.initiate.q[k]), grad(self.eta_diff[j])) * dx - \
                                         self.initiate.stoIJK[i][j][k] * self.initiate.dt * self.initiate.H[i] * \
                                         grad(self.initiate.q[k])[0] * self.ut_mean[2 * j] * dx - \
                                         self.initiate.stoIJK[i][j][k] * self.initiate.dt * self.initiate.H[i] * \
                                         grad(self.initiate.q[k])[1] * self.ut_mean[2 * j + 1] * dx

    def add_convection(self):
        """ add convection term to tentative water velocity F_u_tent weak form. """
        for k in range(self.initiate.n_modes):
            for j in range(self.initiate.n_modes):
                for i in range(self.initiate.n_modes):
                    if abs(self.initiate.stoIJK[i][j][k]) > DOLFIN_EPS:
                        self.F_u_tent += self.initiate.stoIJK[i][j][k] * inner(
                            self.initiate.v[2 * k], self.u_bash[2 * i] * grad(self.u_mean[2 * j])[0] +
                            self.u_bash[2 * i + 1] * grad(self.u_mean[2 * j])[1]) * dx + \
                                         self.initiate.stoIJK[i][j][k] * inner(
                            self.initiate.v[2 * k + 1], self.u_bash[2 * i] * grad(self.u_mean[2 * j + 1])[0] +
                            self.u_bash[2 * i + 1] * grad(self.u_mean[2 * j + 1])[1]) * dx

    def add_viscosity(self):
        """ add viscosity term to tentative water velocity F_u_tent weak form. """
        for k in range(self.initiate.n_modes):
            for j in range(self.initiate.n_modes):
                for i in range(self.initiate.n_modes):
                    if abs(self.initiate.stoIJK[i][j][k]) > DOLFIN_EPS:
                        self.F_u_tent += (self.initiate.stoIJK[i][j][k] * self.initiate.nu_expression_object_list[i] *
                                          inner(grad(self.initiate.v[2 * k]), grad(self.u_mean[2 * j])) * dx) + \
                                         (self.initiate.stoIJK[i][j][k] * self.initiate.nu_expression_object_list[i] *
                                          inner(grad(self.initiate.v[2 * k + 1]), grad(self.u_mean[2 * j + 1])) * dx)

    def add_bottom_stress(self):
        """ add bottom friction source term to tentative water velocity F_u_tent weak form. """
        for k in range(self.initiate.n_modes):
            for j in range(self.initiate.n_modes):
                for i in range(self.initiate.n_modes):
                    if abs(self.initiate.stoIJK[i][j][k]) > DOLFIN_EPS:
                        self.F_u_tent += (self.initiate.stoIJK[i][j][k] * (self.initiate.u0[0] ** 2 +
                                                                           self.initiate.u0[1] ** 2) ** 0.5 /
                                          self.initiate.H[0] * self.initiate.bottomDrag_expression_object_list[i] *
                                          inner(self.u_mean[2 * j], self.initiate.v[2 * k]) * dx) + \
                                         (self.initiate.stoIJK[i][j][k] * (self.initiate.u0[0] ** 2 +
                                                                           self.initiate.u0[1] ** 2) ** 0.5 /
                                          self.initiate.H[0] * self.initiate.bottomDrag_expression_object_list[i] *
                                          inner(self.u_mean[2 * j + 1], self.initiate.v[2 * k + 1]) * dx)

    def add_wind_stress(self):
        """ add wind stress source term to tentative water velocity F_u_tent weak form. """
        norm_wind = (self.initiate.wind_para_x ** 2 + self.initiate.wind_para_y ** 2) ** 0.5
        for k in range(self.initiate.n_modes):
            self.F_u_tent -= (self.initiate.v[2 * k] * self.initiate.rho_air / self.initiate.rho_water * norm_wind *
                              self.initiate.wind_para_x / self.initiate.H[0] *
                              self.initiate.windDrag_expression_object_list[k] * (0.75 + 0.067 * self.norm_wind)) * dx
            self.F_u_tent -= (self.initiate.v[2 * k + 1] * self.initiate.rho_air / self.initiate.rho_water * norm_wind *
                              self.initiate.wind_para_y / self.initiate.H[0] *
                              self.initiate.windDrag_expression_object_list[k] * (0.75 + 0.067 * self.norm_wind)) * dx

    def add_atmospheric_pressure(self):
        """ add atmospheric pressure source term to tentative water velocity F_u_tent weak form. """
        for k in range(self.initiate.n_modes):
            self.F_u_tent += inner(self.initiate.v[2 * k],
                                   grad(self.initiate.pressure / self.inputs.rho_water)[0]) * dx + \
                             inner(self.initiate.v[2 * k + 1],
                                   grad(self.initiate.pressure / self.inputs.rho_water)[1]) * dx

    def add_su_pg(self):
        # TODO : add stochastic version SUPG
        pass

    def add_crosswind(self):
        # TODO: add stochastic version crosswind diffusion
        pass
