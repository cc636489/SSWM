""" class of reading input. """

from input_generalized import *


class ModelReadInput:
    """
    Model Read Input Class.

    Parameters
    ----------
    input_dir : str
        The directory of input files.

    output_dir : str
        The directory of output files.

    output_str : str
        name string of output files.

    n_sample : int
        number of sampling points in random space due to uncertain inputs.

    test_node_x : float
        x coordinate of a physical point in the domain, selected for statistical analysis.

    test_node_y : float
        y coordinate of a physical point in the domain, selected for statistical analysis.

    dist_name : str
        the name of distribution, either "uniform" or "Gaussian".

    sto_ploy_deg : int
        the degree of stochastic polynomial chaos basis.

    sto_poly_dim : int
        the dimension of stochastic polynomial chaos basis.

    coefficient : list
        it represents the lower and upper bound of "uniform" distribution.

    domain : dict
        mesh type and mesh file.

    sto_viscosity : str
        uncertain expression of viscosity. "2*x*q0 + 3*y*q1"

    sto_bottomDrag : str
        uncertain expression of bottom drag(friction) coefficient.

    sto_windDrag : str
        uncertain expression of wind drag coefficient.

    tidal_amplitude : float
        M2 tides amplitude

    tidal_period : float
        M2 tides period

    start_time : float
        model run start time

    end_time : float
        model run end time

    time_step : float
        delta t.

    theta : float
        explicit implicit parameter, 0:explicit, 1:implicit

    bathymetry : dict
        types of bathymetry input, either builtin or import bath file

    initial_u : dict
        types of initial input

    boundary_u : dict
        types of initial input

    """

    def __init__(self):

        self.input_dir = input_dir
        self.output_dir = output_dir
        self.output_str = output_str
        self.wind_file = wind_file
        self.bath_file = bath_file
        self.mesh_file = mesh_file

        self.n_sample = n_sample

        self.test_node_x = test_node_x
        self.test_node_y = test_node_y

        self.dist_name = dist_name
        self.sto_poly_deg = sto_poly_deg
        self.sto_poly_dim = sto_poly_dim
        self.coefficient = coefficient

        self.domain = domain

        self.sto_viscosity = sto_viscosity
        self.sto_bottomDrag = sto_bottomDrag
        self.sto_windDrag = sto_windDrag

        self.include_viscosity = include_viscosity
        self.include_convection = include_convection
        self.linear_divergence = linear_divergence
        self.include_les = include_les
        self.include_wind_stress = include_wind_stress
        self.include_const_wind = include_const_wind
        self.include_bottom_stress = include_bottom_stress
        self.include_atmospheric_pressure = include_atmospheric_pressure
        self.include_supg = include_supg
        self.include_crosswind = include_crosswind
        self.include_auxiliary_viscosity = include_auxiliary_viscosity
        self.include_interior_penalty = include_interior_penalty
        self.les_parameters = les_parameters
        self.DEBUG_mode = DEBUG_mode  # not fully implemented yet
        self.USE_pvd = USE_pvd
        self.USE_HDF5 = USE_HDF5
        self.USE_iterative = USE_iterative

        self.wind_x = wind_x
        self.wind_y = wind_y

        self.tidal_amplitude = tidal_amplitude
        self.tidal_period = tidal_period
        self.start_time = start_time
        self.end_time = end_time
        self.time_step = time_step
        self.theta = theta

        self.wind_dt = wind_dt
        self.R_earth = R_earth
        self.DEG2RAD = DEG2RAD
        self.MPE_deg = MPE_deg
        self.OMEGA = OMEGA
        self.P_background = P_background
        self.bla_dj = bla_dj
        self.rho_air = rho_air
        self.rho_water = rho_water
        self.Gravity = Gravity
        self.one2ten = one2ten
        self.wind_scheme = "powell"
        self.sigma = 9.0 #ike
        self.sigma = 4.0 #harvey

        self.bathymetry = bathymetry

        self.initial_u = initial_u
        self.initial_eta = initial_eta

        self.bc_file = bc_file
        self.boundary_u = boundary_u
        self.boundary_eta = boundary_eta
