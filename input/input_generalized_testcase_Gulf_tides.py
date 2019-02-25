##################################################################
# Input parameter
##################################################################
# output file name string
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
output_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/**_regression_test/"
output_str = "gulf_tides_"
bath_file = "bathymetry_small_gulf.nc"
mesh_file = "fort.14_small_gulf.xml"
wind_file = "fort.22_ike_short_short"
boundary_file = "fort.14_small_gulf_facet_region.xml"

# stochastic input
n_sample = 1
test_node_x = [1.2643567178100001e+06]
test_node_y = [3.3574974964899998e+06]

# stochastic basis
dist_name = "uniform"  # "uniform"  or "gaussian"   ==> only enable "uniform" mode at present.
sto_poly_deg = 0  # polynomial chaos order is 2.
# the dimension and coefficient should be paired and dim!=0  and coefficient!=Null
sto_poly_dim = 2  # use "q0","q1","q2", ....
coefficient = [1, 2, 1, 2]  # lower1/upper1--lower2/upper2--...

# horizontal domain setting
# first way: simple domain built in
# domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
# second way: complex domain import: domain = {"importfile": "inlet.xml"}
domain = {"importfile": input_dir + mesh_file}

# stochastic coefficient # if contains sin(), cos(), should use sympy sin and sympy cos!!!!!!!!!!
sto_viscosity = "1e-6"
sto_bottomDrag = "0.003"
sto_windDrag = "0.0"

# terms control
include_viscosity = True
include_convection = True
linear_divergence = False
include_les = True
include_wind_stress = False
include_bottom_stress = True
include_atmospheric_pressure = False
include_supg = False
include_crosswind = False
les_parameters = {'smagorinsky_coefficient': 1.0}
DEBUG_mode = False
USE_pvd = False
USE_HDF5 = False
USE_iterative = False

# time parameter setting
tidal_amplitude = 0.3
tidal_period = 12.41666*60*60
start_time = 0.0
end_time = 223500.0
time_step = 447.0
theta = 1.0

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
# bathymetry = {"flat": 20.}
# bathymetry = {"vary": "2*x+3*y"}
# bathymetry = {"class": ['type1',400,600,'5.0','-3*((x-500)/100)**4 + 6*((x-500)/100)**2 + 2']}
# bathymetry = {"class": ['type2', 2150, '-14.0/2150.0*x + 19.0', '5.0']}
bathymetry = {"class": ['type3', input_dir + bath_file]}

# initial condition
initial_u = {"flat": (0., 0.)}
initial_eta = {"flat": 0.0}
# initial_u ={"vary": ('2*x+3*y', '3*y')}
# initial_eta = {"vary": "sp.cos(pi/100.0*x) "}

# boundary condition
# first way:
bc_file = input_dir + boundary_file
# boundary_u = {1: "no_slip", 2: "no_slip", 3: "no_slip"}
boundary_u = {}
boundary_eta = {3: "M2 special"}
# boundary_eta = {}
# second way:
# bc_file = None
# boundary_u = {1: "free_slip_in_y", 2: "free_slip_in_y", 3: "free_slip_in_x", 4: "free_slip_in_x"}
# include sin cos ==> should add sympy symbol!!!
# boundary_eta = {}
