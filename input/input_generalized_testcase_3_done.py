##################################################################
# Input parameter
##################################################################
# output file name string
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
output_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test/"
output_str = "test3_generalized_"
bath_file = "hump_adcirc_sswm.nc"
mesh_file = "hump_adcirc_sswm.xml"
wind_file = " "
boundary_file = "hump_adcirc_sswm_facet_region.xml"

# stochastic input
n_sample = 1
test_node_x = [500.]
test_node_y = [100.]

# stochastic basis
dist_name = "uniform"  # "uniform"  or "gaussian"   ==> only enable "uniform" mode at present.
sto_poly_deg = 0  # polynomial chaos order is 2.
# the dimension and coefficient should be paired and dim!=0  and coefficient!=Null
sto_poly_dim = 2  # use "q0","q1","q2", ....
coefficient = [0.8, 1.2, 0.9, 1.1]  # lower1/upper1--lower2/upper2--...

# horizontal domain setting
# first way: simple domain built in
# domain = {"rectangle": (0., 0., 1000., 200., 100, 20)}
# second way: complex domain import: domain = {"importfile": "inlet.xml"}
domain = {"importfile": input_dir + mesh_file}

# stochastic coefficient # if contains sin(), cos(), should use sympy sin and sympy cos!!!!!!!!!!
sto_viscosity = "0.0"
sto_bottomDrag = "0.0"
sto_windDrag = "0.0"

# terms control
include_viscosity = True
include_convection = True
linear_divergence = False
include_les = False
include_wind_stress = False
include_bottom_stress = False
include_atmospheric_pressure = False
include_supg = False
include_crosswind = False
les_parameters = {'smagorinsky_coefficient': 0.13}
DEBUG_mode = False
USE_pvd = False
USE_HDF5 = False
USE_iterative = False

# time parameter setting
tidal_amplitude = 0.75
tidal_period = 12.41666*60*60
start_time = 0.0
end_time = 300.0
time_step = 1.0
theta = 0.5

# equation parameter setting
wind_dt = 3600 * 6
R_earth = 6378206.4
DEG2RAD = 3.1415926 / 180.
MPE_deg = R_earth * DEG2RAD
OMEGA = 2 * 3.1415926 / 86164.2
P_background = 1013.0
bla_dj = 0.9
rho_air = 1.15
rho_water = 1000
Gravity = 9.80665
one2ten = 0.8928


# bathymetry setting
#bathymetry = {"flat": 20.}
#bathymetry = {"vary": "2*x+3*y"}
# bathymetry = {"class": ['type1', 400, 600, '5.0', '-3.0 * ((x-500.0)/100.0)**4 + 6.0 * ((x-500.0)/100.0)**2 + 2.0']}
# bathymetry = {"class": ['type2', 2150, '-14.0/2150.0*x + 19.0', '5.0']}
bathymetry = {"class": ['type3', input_dir+bath_file]}


# initial condition
initial_u = {"flat": (0., 0.)}
initial_eta = {"flat": 0.0}
#initial_u ={"vary": ('2*x+3*y', '3*y')}
#initial_eta = {"vary": "sp.cos(pi/100.0*x) "}

# boundary condition
# first way:
bc_file = input_dir + boundary_file
# boundary_u = {1: "free_slip_in_y", 3: "free_slip_in_x", 4: "free_slip_in_x"}
boundary_u = {}
boundary_eta = {3: "0.1 * sp.sin(pi/20.0 * t)"}
# second way:
# bc_file = None
# boundary_u = {1: "free_slip_in_y", 2: "free_slip_in_y", 3: "free_slip_in_x", 4: "free_slip_in_x"}
# include sin cos ==> should add sympy symbol!!!
# boundary_eta = {}
