##################################################################
# Input parameter
##################################################################
# output file name string
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-NS-github/input/"
outputdir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/test2_res/"
outputstr = "##_test2_generalized_sloth_freeslip_success_"
bathyfile = "bathymetry_small_gulf.nc"
meshfile = "inlet_coarse_50.xml"
windfile = "fort.22_ike_short_short"
boundaryfile = "inlet_coarse_50_facet_region.xml"

# stochastic input
nsample = 1
xytickstep = 4
xyticknumber = 6
scatterplotsize = 5
# only implement with only one testnode at one time.
# test_node_x = [1.2643567178100001e+06]
# test_node_y = [3.3574974964899998e+06]

testnodeX = [25.]
testnodeY = [25.]

# zmax = 1.3
# zmin = -8.3

# stochastic basis
distname = "uniform"  # "uniform"  or "gaussian"   ==> only enable "uniform" mode at present.
sto_poly_deg = 0 # polynomial chaos order is 2.
# the dimension and coefficient should be paired and dim!=0  and coefficient!=Null
sto_poly_dim = 2 # use "q0","q1","q2", ....
coefficient = [-1.2, 1.2, -2, 2] # lower1/upper1--lower2/upper2--...

# horizontal domain setting
# first way: simple domain built in
domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
# second way: complex domain import: domain = {"importfile": "inlet.xml"}
# domain = {"importfile": input_dir+mesh_file}

# stochastic coefficient # if contains sin(), cos(), should use sympy sin and sympy cos!!!!!!!!!!
sto_viscosity = "0.0"
sto_bottomDrag = "0.0"
sto_windDrag = "0.0"

# terms control
include_viscosity = True
include_advection = True
linear_divergence = False
include_les = False
include_wind = False
les_parameters = {'smagorinsky_coefficient': 1.0}
DEBUG_mode = False
USEpvd = True


# time parameter setting
tidal_amplitude = 1.5  # galveston bay area normal tides.
tidal_period = 12*60*60
startime = 0.0
endtime = 50.0
timestep = 0.5
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
bathymetry = {"flat": 20.}
# bathymetry = {"vary": "2*x+3*y"}
# bathymetry = {"class": ['type1',400,600,'5.0','-3*((x-500)/100)**4 + 6*((x-500)/100)**2 + 2']}
# bathymetry = {"class": ['type2', 2150, '-14.0/2150.0*x + 19.0', '5.0']}
# bathymetry = {"class": ['type3', input_dir+bath_file]}

# initial condition
initial_u = {"flat": (0., 0.)}
# initial_eta = {"flat": 0.0}
# initial_u ={"vary": ('2*x+3*y', '3*y')}
initial_eta = {"vary": "0.1*sp.cos(pi/100.0*x)"}

# boundary condition
# first way:
# bc_file = input_dir + boundary_file
# boundary_u = {1: "noslip", 2: "noslip", 3: "noslip"}
boundary_u = {"Left": "freeslipyy", "Right": "freeslipyy", "Up": "freeslipxx", "Down": "freeslipxx"}#{1: "freeslipyy", 2: "freeslipxx"}# # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
# boundary_eta = {3: "inlet special"} # "float/str"
boundary_eta = {}
# second way:
bcfile = None
# boundary_u = {"Left": "freeslipyy", "Right": "freeslipyy", "Up": "freeslipxx", "Down": "freeslipxx"}# include sin cos ==> should add sympy symbol!!! # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
# boundary_eta = {}# "float/str"
