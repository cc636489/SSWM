

from input.input_generalized import *


class ModelReadInput:

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
        self.include_bottom_stress = include_bottom_stress
        self.include_atmospheric_pressure = include_atmospheric_pressure
        self.include_supg = include_supg
        self.include_crosswind = include_crosswind
        self.les_parameters = les_parameters
        self.DEBUG_mode = DEBUG_mode  # not implemented yet
        self.USE_pvd = USE_pvd

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

        self.bathymetry = bathymetry

        self.initial_u = initial_u
        self.initial_eta = initial_eta

        self.bc_file = bc_file
        self.boundary_u = boundary_u
        self.boundary_eta = boundary_eta
