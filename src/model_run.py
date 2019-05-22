

from fenics import lhs, rhs, assemble, LUSolver, KrylovSolver, File, \
                   XDMFFile, HDF5File, dx, solve, Constant, project, as_vector, set_log_level
import numpy as np
set_log_level(30)

class ModelRun:
    """
    class to solve the formulated weak form at each time step.

    Parameters:
    ----------
    inputs : object
        object of class ModelReadInput from model_read_input.py.

    initiate : object
        object of class ModelInitiate from model_initiate.py.

    formulate : object
        object of class (Det)(Sto)ModelFormulate from model_formulate.py.

    a_u_tent : object
        left hand side of F_u_tent

    L_u_tent : object
        right hand side of F_u_tent

    a_p_corr : object
        left hand side of F_p_corr

    L_p_corr : object
        right hand side of F_p_corr

    a_u_corr : object
        left hand side of F_u_corr

    L_u_corr : object
        right hand side of F_u_corr

    test_nodes : nested list
        test node x-y coordinates constructing in a list

    bin_random_u1: ndarray
        array to save u-water velocity at each stochastic samples at each time step at a prescribed phycial point.

    bin_random_v1: ndarray
        array to save v-water velocity at each stochastic samples at each time step at a prescribed phycial point.

    bin_random_eta1: ndarray
        array to save surface elevation at each stochastic samples at each time step at a prescribed phycial point.

    total_mass : list
        array to save total volume mass of the entire domain at each time step.

    time_stamp : list
        array to save time stamp at each time step.

    """

    def __init__(self, inputs, initiate, formulate):

        self.inputs = inputs
        self.initiate = initiate
        self.formulate = formulate

        self.a_u_tent = None
        self.L_u_tent = None
        self.a_p_corr = None
        self.L_p_corr = None
        self.a_u_corr = None
        self.L_u_corr = None

        self.A_p_corr = None
        self.a_p_corr_solver = None

        self.A_u_corr = None
        self.a_u_corr_solver = None

        self.u_file = []
        self.eta_file = []
        self.wind_xy_output_file = []
        self.wind_drag_output_file = []
        self.u_used_for_read_back = []
        self.eta_used_for_read_back = []
        self.wind_xy_used_for_read_back = []
        self.wind_drag_used_for_read_back = []

        self.test_nodes = None
        self.sample_x = None
        self.sample_y = None
        self.bin_random_u1 = None
        self.bin_random_v1 = None
        self.bin_random_eta1 = None

        self.total_mass = []
        self.time_stamp = []

        self.time_step_count = 0

    def run_prep(self):
        """ prepare left and right hand side weak form and do necessary assembling before simulation run. """
        self.a_u_tent = lhs(self.formulate.F_u_tent)
        self.L_u_tent = rhs(self.formulate.F_u_tent)
        self.a_p_corr = lhs(self.formulate.F_p_corr)
        self.L_p_corr = rhs(self.formulate.F_p_corr)
        self.a_u_corr = lhs(self.formulate.F_u_corr)
        self.L_u_corr = rhs(self.formulate.F_u_corr)

        if self.inputs.linear_divergence:
            self.A_p_corr = assemble(self.a_p_corr)
            for bc in self.initiate.eta_bc_object_list:
                bc.apply(self.A_p_corr)
            self.a_p_corr_solver = LUSolver(self.A_p_corr)
            self.a_p_corr_solver.parameters["reuse_factorization"] = True
        else:
            self.a_p_corr_solver = KrylovSolver("gmres", "ilu")
            self.a_p_corr_solver.parameters['absolute_tolerance'] = 1E-4
            self.a_p_corr_solver.parameters['relative_tolerance'] = 1E-3
            self.a_p_corr_solver.parameters['maximum_iterations'] = 1000
            self.a_p_corr_solver.parameters['monitor_convergence'] = True

        self.A_u_corr = assemble(self.a_u_corr)
        for bc in self.initiate.u_bc_object_list:
            bc.apply(self.A_u_corr)
        self.a_u_corr_solver = LUSolver(self.A_u_corr)
        self.a_u_corr_solver.parameters["reuse_factorization"] = True

        self._run_create_file()
        self._run_write_file()
        self._run_initialize_bins()  # include create_bins and zeros the initial time step.
        self._update_time()

    def running(self):
        """ iterate each time step to solve formulated weak forms. """
        while float(self.initiate.t - self.initiate.finish_time) < -1e-3:

            print "starting step:", self.time_step_count, "  time is: ", float(self.initiate.t)

            self._update_boundary()

            if self.inputs.include_wind_stress or self.inputs.include_atmospheric_pressure:
                print "Compute wind stress, pressure and wind drag coefficient."
                self._update_wind()

            if self.inputs.include_const_wind:
                pass

            if self.inputs.include_les:
                print "Compute eddy viscosity."
                self._update_les()

            print "Compute tentative velocity."
            self._update_u_tent()

            print "Compute pressure correction."
            self._update_eta_corr()

            print "Compute velocity update."
            self._update_u_corr()

            self.initiate.u00.assign(self.initiate.u0)
            self.initiate.u0.assign(self.initiate.u1)
            self.initiate.eta0.assign(self.initiate.eta1)
            print("===========================================")

            self._run_write_file()
            self._run_write_bins()

            #import pdb;pdb.set_trace()
            self._update_time()

    def run_final(self):
        """ save results to files. """
        self._run_save_bin()
        print "total time_step is ", self.time_step_count

    def _run_create_file(self):
        """ open files to write results. """
        for mode in range(self.initiate.n_modes):
            if self.inputs.USE_pvd:
                self.u_file.append(File(self.inputs.output_dir + "u_" + self.inputs.output_str +
                                        "{:02d}".format(mode) + ".pvd"))
                self.eta_file.append(File(self.inputs.output_dir + "eta_" + self.inputs.output_str +
                                          "{:02d}".format(mode) + ".pvd"))
            else:
                self.u_file.append(XDMFFile(self.inputs.output_dir + "u_" + self.inputs.output_str +
                                            "{:02d}".format(mode) + ".xdmf"))
                self.eta_file.append(XDMFFile(self.inputs.output_dir + "eta_" + self.inputs.output_str +
                                              "{:02d}".format(mode) + ".xdmf"))

            if self.inputs.USE_HDF5:
                self.u_used_for_read_back.append(HDF5File(self.initiate.mesh.mpi_comm(), self.inputs.output_dir +
                                                          "u_used_for_read_back_" + self.inputs.output_str +
                                                          "{:02d}".format(mode) + ".h5", "w"))
                self.eta_used_for_read_back.append(HDF5File(self.initiate.mesh.mpi_comm(), self.inputs.output_dir +
                                                            "eta_used_for_read_back_" + self.inputs.output_str +
                                                            "{:02d}".format(mode) + ".h5", "w"))

        if self.inputs.include_wind_stress:
            if self.inputs.USE_pvd:
                self.wind_xy_output_file.append(File(self.inputs.output_dir + "wind_vel_" +
                                                     self.inputs.output_str + "0.pvd"))
                self.wind_drag_output_file.append(File(self.inputs.output_dir + "wind_drag_" +
                                                     self.inputs.output_str + "0.pvd"))
            else:
                self.wind_xy_output_file.append(XDMFFile(self.inputs.output_dir + "wind_vel_" +
                                                         self.inputs.output_str + "0.xdmf"))
                self.wind_drag_output_file.append(XDMFFile(self.inputs.output_dir + "wind_drag_" +
                                                         self.inputs.output_str + "0.xdmf"))
            if self.inputs.USE_HDF5:
                self.wind_xy_used_for_read_back.append(HDF5File(self.initiate.mesh.mpi_comm(), self.inputs.output_dir +
                                                                "wind_used_for_read_back_" + self.inputs.output_str +
                                                                "0.h5", "w"))
                self.wind_drag_used_for_read_back.append(HDF5File(self.initiate.mesh.mpi_comm(), self.inputs.output_dir +
                                                                "wind_drag_used_for_read_back_" + self.inputs.output_str +
                                                                "0.h5", "w"))

    def _run_write_file(self):
        """ write to files. """
        if self.initiate.n_modes == 1:
            if self.inputs.USE_pvd:
                self.u_file[0] << (self.initiate.u1, float(self.initiate.t))
                self.eta_file[0] << (self.initiate.eta1, float(self.initiate.t))
            else:
                self.u_file[0].write(self.initiate.u1, float(self.initiate.t))
                self.eta_file[0].write(self.initiate.eta1, float(self.initiate.t))
            if self.inputs.USE_HDF5:
                self.u_used_for_read_back[0].write(self.initiate.u1, "WaterVelocity", float(self.initiate.t))
                self.eta_used_for_read_back[0].write(self.initiate.eta1, "SurfaceElevation", float(self.initiate.t))
        else:
            for mode in range(self.initiate.n_modes):
                tmp1 = self.initiate.u1.split(deepcopy=True)
                tmp2 = self.initiate.eta1.split(deepcopy=True)
                if self.inputs.USE_pvd:
                    self.u_file[mode] << (tmp1[mode], float(self.initiate.t))
                    self.eta_file[mode] << (tmp2[mode], float(self.initiate.t))
                else:
                    self.u_file[mode].write(tmp1[mode], float(self.initiate.t))
                    self.eta_file[mode].write(tmp2[mode], float(self.initiate.t))
                if self.inputs.USE_HDF5:
                    self.u_used_for_read_back[mode].write(tmp1[mode], "WaterVelocity", float(self.initiate.t))
                    self.eta_used_for_read_back[mode].write(tmp2[mode], "SurfaceElevation", float(self.initiate.t))

        if self.inputs.include_wind_stress:
            self.initiate.wind_vector.assign(project(as_vector([self.initiate.wind_para_x, self.initiate.wind_para_y]),
                                                     self.initiate.C))
            self.initiate.wind_drag.assign(project(self.initiate.wdrag, self.initiate.B))
            if self.inputs.USE_pvd:
                self.wind_xy_output_file[0] << (self.initiate.wind_vector, float(self.initiate.t))
                self.wind_drag_output_file[0] << (self.initiate.wind_drag, float(self.initiate.t))
            else:
                self.wind_xy_output_file[0].write(self.initiate.wind_vector, float(self.initiate.t))
                self.wind_drag_output_file[0].write(self.initiate.wind_drag, float(self.initiate.t))
            if self.inputs.USE_HDF5:
                self.wind_xy_used_for_read_back[0].write(self.initiate.wind_vector, "WindSpeed", float(self.initiate.t))
                self.wind_drag_used_for_read_back[0].write(self.initiate.wind_drag, "WindDrag", float(self.initiate.t))

    def _run_initialize_bins(self):
        """ initialize arrays to save stochastic results. """
        if self.initiate.n_modes == 1:
            self.total_mass.append(assemble(self.initiate.H * dx))
        self.time_stamp.append(float(self.initiate.t))

    def _run_write_bins(self):
        """ write to array of stochastic results. """
        if self.initiate.n_modes == 1:
            self.total_mass.append(assemble(self.initiate.H * dx))
        self.time_stamp.append(float(self.initiate.t))

    def _run_save_bin(self):
        """ save to file of the arrays. """
        np.save(self.inputs.output_dir + "total_mass_at_every_time_step", self.total_mass)
        np.save(self.inputs.output_dir + "time_stamp_at_every_time_step", self.time_stamp)

    def _update_time(self):
        """ update current time. """
        self.initiate.t = Constant(self.initiate.t + self.initiate.dt)
        self.time_step_count += 1

    def _update_boundary(self):
        """ update boundary condition if it's time dependent. """
        if self.initiate.uTimeDependent:
            self.initiate.u_list_expression.t = self.initiate.t
        if self.initiate.etaTimeDependent:
            self.initiate.eta_list_expression.t = self.initiate.t

    def _update_wind(self):
        """ update wind field if it's time dependent. """
        self.initiate.wind.get_prepared(current_time=float(self.initiate.t))
        wind_para_x_list = []
        wind_para_y_list = []
        pressure_list = []
        if self.inputs.wind_scheme == "powell":
            wdrag_list = []
        for ts in range(len(self.initiate.x_deg)):
            wind_x, wind_y, press = self.initiate.wind.get_wind_x_wind_y(
                    x_coord_deg=self.initiate.x_deg[ts], y_coord_deg=self.initiate.y_deg[ts])
            wind_para_x_list.append(wind_x)
            wind_para_y_list.append(wind_y)
            pressure_list.append(press)
            if self.inputs.wind_scheme == "powell":
                drag = self.initiate.wind.get_wind_drag(x_coord_deg=self.initiate.x_deg[ts], 
                                                        y_coord_deg=self.initiate.y_deg[ts])
                wdrag_list.append(drag)
        for sm in range(len(self.initiate.x_coord)):
            self.initiate.wind_para_x.update(self.initiate.x_coord[sm], self.initiate.y_coord[sm], wind_para_x_list[sm])
            self.initiate.wind_para_y.update(self.initiate.x_coord[sm], self.initiate.y_coord[sm], wind_para_y_list[sm])
            self.initiate.pressure.update(self.initiate.x_coord[sm], self.initiate.y_coord[sm], pressure_list[sm])
            if self.inputs.wind_scheme == "powell":
                self.initiate.wdrag.update(self.initiate.x_coord[sm], self.initiate.y_coord[sm], wdrag_list[sm])

    def _update_les(self):
        """ update eddy viscosity field. """
        if self.initiate.n_modes == 1:
            self.initiate.les[0].u = self.initiate.u0
            self.initiate.les[0].solve()
        else:
            for k in range(self.initiate.n_modes):
                self.initiate.les[k].u = self.initiate.u0.split()[k]
                self.initiate.les[k].solve()

    def _update_u_tent(self):
        """ update form for tentative velcoity. """
        aa_u_tent = assemble(self.a_u_tent)
        b = assemble(self.L_u_tent)
        for bc in self.initiate.u_bc_object_list:
            bc.apply(aa_u_tent, b)
        if self.inputs.USE_iterative:
            solve(aa_u_tent, self.initiate.ut.vector(), b, "gmres", "default")
        else:
            solve(aa_u_tent, self.initiate.ut.vector(), b)

    def _update_eta_corr(self):
        """ update form for the (n+1)th time step pressure. """
        b = assemble(self.L_p_corr)
        for bc in self.initiate.eta_bc_object_list:
            bc.apply(b)
        if self.inputs.linear_divergence:
            self.a_p_corr_solver.solve(self.initiate.eta1.vector(), b)
        else:
            aa_p_corr = assemble(self.a_p_corr)
            for bc in self.initiate.eta_bc_object_list:
                bc.apply(aa_p_corr)
            solve(aa_p_corr, self.initiate.eta1.vector(), b)

    def _update_u_corr(self):
        """ update form for the (n+1)th time step velocity. """
        b = assemble(self.L_u_corr)
        for bc in self.initiate.u_bc_object_list:
            bc.apply(b)
        self.a_u_corr_solver.solve(self.initiate.u1.vector(), b)
