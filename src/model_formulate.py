

from fenics import inner, grad, DOLFIN_EPS, dx, ds, dS, \
    CellSize, div, dot, conditional, gt, FacetNormal, FacetArea, jump


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

        if self.inputs.include_wind_stress:
            self.norm_wind = (self.initiate.wind_para_x ** 2 + self.initiate.wind_para_y ** 2) ** 0.5

        self.F_u_tent = 0
        self.F_p_corr = 0
        self.F_u_corr = 0

        self.h = CellSize(self.initiate.mesh)
        # self.tau1 = 0.5*h*pow(4.0*nu/h+2.0*u_norm, -1.0)
        # self.tau1 = pow((1./0.5/self.initiate.dt)**2 + (2.*self.u0_norm/self.h)**2 +
        # 9.0*(4.*self.initiate.nu_expression_object_list[0]/self.h/self.h)**2, -1.0/2)
        self.tau1 = pow(1 / 0.5 / self.initiate.dt + 2. * self.u0_norm / self.h +
                        4. * self.initiate.nu_expression_object_list[0] / self.h / self.h, -1.0)
        self.tau2_u = self.tau1
        self.tau2_v = self.tau1

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
        self.F_u_tent += self.initiate.nu_expression_object_list[0] * inner(grad(self.initiate.v),
                                                                            grad(self.u_mean)) * dx

    def add_bottom_stress(self):
        """ add bottom friction source term to tentative water velocity F_u_tent weak form. """
        self.F_u_tent += self.u0_norm / self.initiate.H * self.initiate.bottomDrag_expression_object_list[0] * inner(
            self.u_mean, self.initiate.v) * dx

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
                  self.inputs.rho_water * self.norm_wind * self.initiate.wind_para_x / self.initiate.H *
                  self.initiate.windDrag_expression_object_list[0] * (0.75 + 0.067 * self.norm_wind)) * dx
            r -= (self.tau1 * dot(self.u_bash, grad(self.initiate.v[1])) * self.inputs.rho_air /
                  self.inputs.rho_water * self.norm_wind * self.initiate.wind_para_y / self.initiate.H *
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

        u_bash_parallel_norm = (u_bash_parallel[0] ** 2 + u_bash_parallel[1] ** 2) ** 0.5
        v_bash_parallel_norm = (v_bash_parallel[0] ** 2 + v_bash_parallel[1] ** 2) ** 0.5

        tau_parallel_u = pow(1 / 0.5 / self.initiate.dt + 2. * u_bash_parallel_norm / self.h +
                             4. * self.initiate.nu_expression_object_list[0] / self.h / self.h, -1.0)
        tau_parallel_v = pow(1 / 0.5 / self.initiate.dt + 2. * v_bash_parallel_norm / self.h +
                             4. * self.initiate.nu_expression_object_list[0] / self.h / self.h, -1.0)

        # noticing that, this self.tau2 setting assumes the mesh is not distorted.
        self.tau2_u = conditional(gt(tau_parallel_u - self.tau1, 0.0), tau_parallel_u - self.tau1, 0.0)
        self.tau2_v = conditional(gt(tau_parallel_u - self.tau1, 0.0), tau_parallel_v - self.tau1, 0.0)

        ru = inner((1 / self.initiate.dt) * self.u_diff[0] + self.initiate.g * grad(self.initiate.eta0)[0],
                   self.tau2_u * dot(grad(self.initiate.v[0]), u_bash_parallel)) * dx
        rv = inner((1 / self.initiate.dt) * self.u_diff[1] + self.initiate.g * grad(self.initiate.eta0)[1],
                   self.tau2_v * dot(grad(self.initiate.v[1]), v_bash_parallel)) * dx

        if self.inputs.include_convection:
            ru += inner(dot(grad(self.u_mean[0]), self.u_bash), self.tau2_u *
                        dot(grad(self.initiate.v[0]), u_bash_parallel)) * dx
            rv += inner(dot(grad(self.u_mean[1]), self.u_bash), self.tau2_v *
                        dot(grad(self.initiate.v[1]), v_bash_parallel)) * dx

        if self.inputs.include_viscosity:
            ru += inner(-div(self.initiate.nu_expression_object_list[0] * grad(self.u_mean[0])),
                        self.tau2_u * dot(grad(self.initiate.v[0]), u_bash_parallel)) * dx
            rv += inner(-div(self.initiate.nu_expression_object_list[0] * grad(self.u_mean[1])),
                        self.tau2_v * dot(grad(self.initiate.v[1]), v_bash_parallel)) * dx

        if self.inputs.include_bottom_stress:
            ru += inner(self.u0_norm / self.initiate.H * self.initiate.bottomDrag_expression_object_list[0] *
                        self.u_mean[0], self.tau2_u * dot(grad(self.initiate.v[0]), u_bash_parallel)) * dx
            rv += inner(self.u0_norm / self.initiate.H * self.initiate.bottomDrag_expression_object_list[0] *
                        self.u_mean[1], self.tau2_v * dot(grad(self.initiate.v[1]), v_bash_parallel)) * dx

        if self.inputs.include_wind_stress:
            ru -= (self.tau2_u * dot(grad(self.initiate.v[0]), u_bash_parallel) * self.inputs.rho_air /
                   self.inputs.rho_water * self.norm_wind * self.initiate.wind_para_x / self.initiate.H *
                   self.initiate.windDrag_expression_object_list[0] * (0.75 + 0.067 * self.norm_wind)) * dx
            rv -= (self.tau2_v * dot(grad(self.initiate.v[1]), v_bash_parallel) * self.inputs.rho_air /
                   self.inputs.rho_water * self.norm_wind * self.initiate.wind_para_y / self.initiate.H *
                   self.initiate.windDrag_expression_object_list[0] * (0.75 + 0.067 * self.norm_wind)) * dx

        if self.inputs.include_atmospheric_pressure:
            ru += inner(grad(self.initiate.pressure / self.inputs.rho_water)[0],
                        self.tau2_u * dot(grad(self.initiate.v[0]), u_bash_parallel)) * dx
            rv += inner(grad(self.initiate.pressure / self.inputs.rho_water)[1],
                        self.tau2_v * dot(grad(self.initiate.v[1]), v_bash_parallel)) * dx

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
        self.u0_norm = (self.initiate.u0[0] ** 2 + self.initiate.u0[1] ** 2) ** 0.5

        self.eta_diff = self.initiate.eta - self.initiate.eta0
        self.ut_mean = self.initiate.theta * self.initiate.ut + (1. - self.initiate.theta) * self.initiate.u0

        self.eta1_diff = self.initiate.eta1 - self.initiate.eta0

        if self.inputs.include_wind_stress or self.inputs.include_const_wind:
            self.norm_wind = (self.initiate.wind_para_x ** 2 + self.initiate.wind_para_y ** 2) ** 0.5

        self.F_u_tent = 0
        self.F_p_corr = 0
        self.F_u_corr = 0
        self.R = 0

        self.h = CellSize(self.initiate.mesh)
        self.h_edge = FacetArea(self.initiate.mesh)('+')
        self.n = FacetNormal(self.initiate.mesh)
        # self.tau1 = 0.5*h*pow(4.0*nu/h+2.0*u_norm, -1.0)
        # self.tau1 = pow( (1./0.5/dt)**2 + (2.*u_norm/h)**2 + 9.0*(4.*nu/h/h)**2, -1.0/2)
        # TODO: may consider modify tau2 into max(tau2-tau1, 0), now larger consistent viscosity is fine.
        self.tau1 = pow(1 / 0.5 / self.initiate.dt + 2. * self.u0_norm / self.h +
                        4. * self.initiate.nu_expression_object_list[0] / self.h / self.h, -1.0)
        self.tau2 = self.tau1

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

                        # TODO: try to use various stabilized terms, just for test case 4 and above.currently satisfied.
                        self.F_u_tent += (7 * self.u0_norm * self.h / 2.0 * self.initiate.stoIJK[i][j][k] *
                                          self.initiate.nu_expression_object_list[i] *
                                          inner(grad(self.initiate.v[2 * k]), grad(self.u_mean[2 * j])) * dx) + \
                                         (7 * self.u0_norm * self.h / 2.0 * self.initiate.stoIJK[i][j][k] *
                                          self.initiate.nu_expression_object_list[i] *
                                          inner(grad(self.initiate.v[2 * k + 1]), grad(self.u_mean[2 * j + 1])) * dx)

                        # TODO: Nitsche type method: (u_bash,n) = 0 on boundary, so adding zeros to the system.
                        # self.F_u_tent += ((self.initiate.u0[2 * i] * self.n[0] + self.initiate.u0[2 * i + 1] * self.n[1]) *
                        #                   self.initiate.stoIJK[i][j][k] *
                        #                   inner(self.initiate.v[2 * k], self.u_mean[2 * j]) * ds) + \
                        #                  ((self.initiate.u0[2 * i] * self.n[0] + self.initiate.u0[2 * i + 1] * self.n[1]) *
                        #                   self.initiate.stoIJK[i][j][k] *
                        #                   inner(self.initiate.v[2 * k + 1], self.u_mean[2 * j + 1]) * ds)

                        # TODO: Interior penalty method from OpenTidalFarm, if smooth, should be zeros:
                        self.F_u_tent += (1.0 * self.initiate.nu_expression_object_list[i] / self.h_edge *
                                          self.initiate.stoIJK[i][j][k] *
                                          inner(jump(grad(self.initiate.v[2 * k]), self.n),
                                                jump(grad(self.u_mean[2 * j]), self.n)) * dS) + \
                                         (1.0 * self.initiate.nu_expression_object_list[i] / self.h_edge *
                                          self.initiate.stoIJK[i][j][k] *
                                          inner(jump(grad(self.initiate.v[2 * k + 1]), self.n),
                                                jump(grad(self.u_mean[2 * j + 1]), self.n)) * dS)

                        # TODO: Interior penalty method from yuxianglin's forwarded code, if smooth, should be zeros:
                        # self.F_u_tent += (1.0 * self.h_edge ** 2.0 *
                        #                   self.initiate.stoIJK[i][j][k] *
                        #                   inner(jump(grad(self.initiate.v[2 * k]), self.n),
                        #                         jump(grad(self.u_mean[2 * j]), self.n)) * dS) + \
                        #                  (1.0 * self.h_edge ** 2.0 *
                        #                   self.initiate.stoIJK[i][j][k] *
                        #                   inner(jump(grad(self.initiate.v[2 * k + 1]), self.n),
                        #                         jump(grad(self.u_mean[2 * j + 1]), self.n)) * dS)

    def add_bottom_stress(self):
        """ add bottom friction source term to tentative water velocity F_u_tent weak form. """
        for k in range(self.initiate.n_modes):
            for j in range(self.initiate.n_modes):
                for i in range(self.initiate.n_modes):
                    if abs(self.initiate.stoIJK[i][j][k]) > DOLFIN_EPS:
                        self.F_u_tent += (self.initiate.stoIJK[i][j][k] * self.u0_norm /
                                          self.initiate.H[0] * self.initiate.bottomDrag_expression_object_list[i] *
                                          inner(self.u_mean[2 * j], self.initiate.v[2 * k]) * dx) + \
                                         (self.initiate.stoIJK[i][j][k] * self.u0_norm /
                                          self.initiate.H[0] * self.initiate.bottomDrag_expression_object_list[i] *
                                          inner(self.u_mean[2 * j + 1], self.initiate.v[2 * k + 1]) * dx)

    def add_wind_stress(self):
        """ add wind stress source term to tentative water velocity F_u_tent weak form. """
        for k in range(self.initiate.n_modes):
            self.F_u_tent -= (self.initiate.v[2 * k] * self.inputs.rho_air / self.inputs.rho_water *
                              self.norm_wind * self.initiate.wind_para_x / self.initiate.H[0] *
                              self.initiate.windDrag_expression_object_list[k] * (0.75 + 0.067 * self.norm_wind)) * dx
            self.F_u_tent -= (self.initiate.v[2 * k + 1] * self.inputs.rho_air / self.inputs.rho_water *
                              self.norm_wind * self.initiate.wind_para_y / self.initiate.H[0] *
                              self.initiate.windDrag_expression_object_list[k] * (0.75 + 0.067 * self.norm_wind)) * dx

    def add_atmospheric_pressure(self):
        """ add atmospheric pressure source term to tentative water velocity F_u_tent weak form. """
        for k in range(self.initiate.n_modes):
            self.F_u_tent += inner(self.initiate.v[2 * k],
                                   grad(self.initiate.pressure / self.inputs.rho_water)[0]) * dx + \
                             inner(self.initiate.v[2 * k + 1],
                                   grad(self.initiate.pressure / self.inputs.rho_water)[1]) * dx

    def add_su_pg(self):
        """ add standard stabilized term to tentative water velocity F_u_tent weak form. """
        for k in range(self.initiate.n_modes):
            for j in range(self.initiate.n_modes):
                for i in range(self.initiate.n_modes):
                    if abs(self.initiate.stoIJK[i][j][k]) > DOLFIN_EPS:
                        # add SUPG transient term.
                        r = self.initiate.stoIJK[i][j][k] * self.tau1 / self.initiate.dt * \
                            inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j])[0] +
                                  self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j])[1],
                                  self.u_diff[2 * k]
                                  ) * dx + \
                            self.initiate.stoIJK[i][j][k] * self.tau1 / self.initiate.dt * \
                            inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j + 1])[0] +
                                  self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j + 1])[1],
                                  self.u_diff[2 * k + 1]
                                  ) * dx

                        # add SUPG gravity term.
                        r += self.initiate.stoIJK[i][j][k] * self.tau1 * self.initiate.g * \
                            inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j])[0] +
                                  self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j])[1],
                                  grad(self.initiate.eta0[k])[0]
                                  ) * dx + \
                            self.initiate.stoIJK[i][j][k] * self.tau1 * self.initiate.g * \
                            inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j + 1])[0] +
                                  self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j + 1])[1],
                                  grad(self.initiate.eta0[k])[1]
                                  ) * dx

                        # add SUPG convection term.
                        if self.inputs.include_convection:
                            r += self.initiate.stoIJK[i][j][k] * self.initiate.stoIJK[i][j][k] * self.tau1 * \
                                 inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j])[0] +
                                       self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j])[1],
                                       self.u_bash[2 * i] * grad(self.u_mean[2 * j])[0] +
                                       self.u_bash[2 * i + 1] * grad(self.u_mean[2 * j])[1]
                                       ) * dx + \
                                 self.initiate.stoIJK[i][j][k] * self.initiate.stoIJK[i][j][k] * self.tau1 * \
                                 inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j + 1])[0] +
                                       self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j + 1])[1],
                                       self.u_bash[2 * i] * grad(self.u_mean[2 * j + 1])[0] +
                                       self.u_bash[2 * i + 1] * grad(self.u_mean[2 * j + 1])[1]
                                       ) * dx

                        # add SUPG viscosity term.
                        if self.inputs.include_viscosity:
                            r += self.initiate.stoIJK[i][j][k] * self.initiate.stoIJK[i][j][k] * self.tau1 * \
                                 inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j])[0] +
                                       self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j])[1],
                                       -div(self.initiate.nu_expression_object_list[i] * grad(self.u_mean[2 * j]))
                                       ) * dx + \
                                 self.initiate.stoIJK[i][j][k] * self.initiate.stoIJK[i][j][k] * self.tau1 * \
                                 inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j + 1])[0] +
                                       self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j + 1])[1],
                                       -div(self.initiate.nu_expression_object_list[i] * grad(self.u_mean[2 * j + 1]))
                                       ) * dx

                        # add SUPG bottom friction forcing term.
                        if self.inputs.include_bottom_stress:
                            r += (self.initiate.stoIJK[i][j][k] * self.initiate.stoIJK[i][j][k] * self.tau1 *
                                  self.u0_norm / self.initiate.H[0] *
                                  self.initiate.bottomDrag_expression_object_list[i] *
                                  inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j])[0] +
                                        self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j])[1],
                                        self.u_mean[2 * j])
                                  ) * dx + \
                                 (self.initiate.stoIJK[i][j][k] * self.initiate.stoIJK[i][j][k] * self.tau1 *
                                  self.u0_norm / self.initiate.H[0] *
                                  self.initiate.bottomDrag_expression_object_list[i] *
                                  inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j + 1])[0] +
                                        self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j + 1])[1],
                                        self.u_mean[2 * j + 1])
                                  ) * dx

                        # add SUPG wind forcing term.
                        if self.inputs.include_wind_stress or self.inputs.include_const_wind:
                            r -= (self.initiate.stoIJK[i][j][k] * self.tau1 * self.inputs.rho_air /
                                  self.inputs.rho_water * self.norm_wind / self.initiate.H[0] *
                                  self.initiate.windDrag_expression_object_list[k] * (0.75 + 0.067 * self.norm_wind) *
                                  inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j])[0] +
                                        self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j])[1],
                                        self.initiate.wind_para_x)) * dx
                            r -= (self.initiate.stoIJK[i][j][k] * self.tau1 * self.inputs.rho_air /
                                  self.inputs.rho_water * self.norm_wind / self.initiate.H[0] *
                                  self.initiate.windDrag_expression_object_list[k] * (0.75 + 0.067 * self.norm_wind) *
                                  inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j + 1])[0] +
                                        self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j + 1])[1],
                                        self.initiate.wind_para_y)) * dx

                        # add SUPG atmospheric pressure term.
                        if self.inputs.include_atmospheric_pressure:
                            r += (self.initiate.stoIJK[i][j][k] * self.tau1 / self.inputs.rho_water *
                                  inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j])[0] +
                                        self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j])[1],
                                        grad(self.initiate.pressure)[0])) * dx + \
                                 (self.initiate.stoIJK[i][j][k] * self.tau1 / self.inputs.rho_water *
                                  inner(self.u_bash[2 * i] * grad(self.initiate.v[2 * j + 1])[0] +
                                        self.u_bash[2 * i + 1] * grad(self.initiate.v[2 * j + 1])[1],
                                        grad(self.initiate.pressure)[1])) * dx

                        self.F_u_tent += r

    def add_crosswind(self):
        """ add cross wind diffusion term to tentative water velocity F_u_tent weak form. """
        u0_mean = self.initiate.theta * self.initiate.u0 + (1. - self.initiate.theta) * self.initiate.u00

        for k in range(self.initiate.n_modes):
            for j in range(self.initiate.n_modes):
                du_norm = inner(grad(self.initiate.u0[2 * j]), grad(self.initiate.u0[2 * j])) ** 0.5
                dv_norm = inner(grad(self.initiate.u0[2 * j + 1]), grad(self.initiate.u0[2 * j + 1])) ** 0.5
                for i in range(self.initiate.n_modes):
                    if abs(self.initiate.stoIJK[i][j][k]) > DOLFIN_EPS:
                        u_bash_parallel = (self.u_bash[2 * i] * grad(u0_mean[2 * j])[0] +
                                           self.u_bash[2 * i + 1] * grad(u0_mean[2 * j])[1]
                                           ) * grad(u0_mean[2 * j]) / du_norm ** 2
                        v_bash_parallel = (self.u_bash[2 * i] * grad(u0_mean[2 * j + 1])[0] +
                                           self.u_bash[2 * i + 1] * grad(u0_mean[2 * j + 1])[1]
                                           ) * grad(u0_mean[2 * j + 1]) / dv_norm ** 2

                        # add crosswind transient term.
                        ru = self.initiate.stoIJK[i][j][k] * self.tau2 / self.initiate.dt * \
                            inner(dot(grad(self.initiate.v[2 * j]), u_bash_parallel),
                                  self.u_diff[2 * k]) * dx
                        rv = self.initiate.stoIJK[i][j][k] * self.tau2 / self.initiate.dt * \
                            inner(dot(grad(self.initiate.v[2 * j + 1]), v_bash_parallel),
                                  self.u_diff[2 * k + 1]) * dx

                        # add crosswind gravity term.
                        ru += self.initiate.stoIJK[i][j][k] * self.tau2 * self.initiate.g * \
                            inner(dot(grad(self.initiate.v[2 * j]), u_bash_parallel),
                                  grad(self.initiate.eta0[k])[0]) * dx
                        rv += self.initiate.stoIJK[i][j][k] * self.tau2 * self.initiate.g * \
                            inner(dot(grad(self.initiate.v[2 * j + 1]), v_bash_parallel),
                                  grad(self.initiate.eta0[k])[1]) * dx

                        # add crosswind convective term.
                        if self.inputs.include_convection:
                            ru += self.initiate.stoIJK[i][j][k] * self.initiate.stoIJK[i][j][k] * self.tau2 * \
                                inner(dot(grad(self.initiate.v[2 * j]), u_bash_parallel),
                                      self.u_bash[2 * i] * grad(self.u_mean[2 * j])[0] +
                                      self.u_bash[2 * i + 1] * grad(self.u_mean[2 * j])[1]) * dx
                            rv += self.initiate.stoIJK[i][j][k] * self.initiate.stoIJK[i][j][k] * self.tau2 * \
                                inner(dot(grad(self.initiate.v[2 * j + 1]), v_bash_parallel),
                                      self.u_bash[2 * i] * grad(self.u_mean[2 * j + 1])[0] +
                                      self.u_bash[2 * i + 1] * grad(self.u_mean[2 * j + 1])[1]) * dx

                        # add crosswind viscosity term.
                        if self.inputs.include_viscosity:
                            ru += self.initiate.stoIJK[i][j][k] * self.initiate.stoIJK[i][j][k] * self.tau2 * \
                                inner(dot(grad(self.initiate.v[2 * j]), u_bash_parallel),
                                      -div(self.initiate.nu_expression_object_list[i] * grad(self.u_mean[2 * j]))
                                      ) * dx
                            rv += self.initiate.stoIJK[i][j][k] * self.initiate.stoIJK[i][j][k] * self.tau2 * \
                                inner(dot(grad(self.initiate.v[2 * j + 1]), v_bash_parallel),
                                      -div(self.initiate.nu_expression_object_list[i] * grad(self.u_mean[2 * j + 1]))
                                      ) * dx

                        # add crosswind bottom friction forcing term.
                        if self.inputs.include_bottom_stress:
                            ru += (self.initiate.stoIJK[i][j][k] * self.initiate.stoIJK[i][j][k] * self.tau2 *
                                   self.u0_norm / self.initiate.H[0] *
                                   self.initiate.bottomDrag_expression_object_list[i] *
                                   inner(dot(grad(self.initiate.v[2 * j]), u_bash_parallel), self.u_mean[2 * j])
                                   ) * dx
                            rv += (self.initiate.stoIJK[i][j][k] * self.initiate.stoIJK[i][j][k] * self.tau2 *
                                   self.u0_norm / self.initiate.H[0] *
                                   self.initiate.bottomDrag_expression_object_list[i] *
                                   inner(dot(grad(self.initiate.v[2 * j + 1]), v_bash_parallel), self.u_mean[2 * j + 1])
                                   ) * dx

                        # add crosswind wind forcing term.
                        if self.inputs.include_wind_stress or self.inputs.include_const_wind:
                            ru -= (self.initiate.stoIJK[i][j][k] * self.tau2 * self.inputs.rho_air /
                                   self.inputs.rho_water * self.norm_wind / self.initiate.H[0] *
                                   self.initiate.windDrag_expression_object_list[k] * (0.75 + 0.067 * self.norm_wind) *
                                   inner(dot(grad(self.initiate.v[2 * j]), u_bash_parallel),
                                         self.initiate.wind_para_x)) * dx
                            rv -= (self.initiate.stoIJK[i][j][k] * self.tau2 * self.inputs.rho_air /
                                   self.inputs.rho_water * self.norm_wind / self.initiate.H[0] *
                                   self.initiate.windDrag_expression_object_list[k] * (0.75 + 0.067 * self.norm_wind) *
                                   inner(dot(grad(self.initiate.v[2 * j + 1]), v_bash_parallel),
                                         self.initiate.wind_para_y)) * dx

                        # add crosswind atmospheric pressure term.
                        if self.inputs.include_atmospheric_pressure:
                            ru += (self.initiate.stoIJK[i][j][k] * self.tau2 / self.inputs.rho_water * inner(
                                        dot(grad(self.initiate.v[2 * j]), u_bash_parallel),
                                        grad(self.initiate.pressure)[0])) * dx
                            rv += (self.initiate.stoIJK[i][j][k] * self.tau2 / self.inputs.rho_water * inner(
                                        dot(grad(self.initiate.v[2 * j + 1]), v_bash_parallel),
                                        grad(self.initiate.pressure)[1])) * dx

                        self.F_u_tent += ru + rv
