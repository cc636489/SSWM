##################################################################
# Input parameter
##################################################################
# output file name string
outputdir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/test4_res/scratch2/"
outputstr = "test_quadratic_bottom_friction_"

# stochastic input
nsample = 10
xytickstep = 4
xyticknumber = 6
scatterplotsize = 5
#only implement with only one testnode at one time.
testnodeX = [2250.]
testnodeY = [1500.]
#zmax = 1.3
#zmin = -8.3

# stochastic basis
distname = "uniform"  # "uniform"  or "gaussian"   ==> only enable "uniform" mode at present.
sto_poly_deg = 0 # polynomial chaos order is 2.
# the dimension and coefficient should be paired and dim!=0  and coefficient!=Null
sto_poly_dim = 2 # use "q0","q1","q2", ....
coefficient = [1, 2, 1, 2] # lower1/upper1--lower2/upper2--...

# horizontal domain setting
# first way: simple domain built in
# domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
# second way: complex domain import: domain = {"importfile": "inlet.xml"}
domain = {"importfile": "inlet_coarse_50.xml"}

# stochastic coefficient # if contains sin(), cos(), should use sympy sin and sympy cos!!!!!!!!!!
sto_viscosity = "1e-6"
sto_bottomDrag = "0.003"
sto_windDrag = "0.0013"

# terms control
include_viscosity = True
include_advection = True
linear_divergence = False
include_les = True
les_parameters = {'smagorinsky_coefficient': 1.0}
DEBUG_mode = False
USEpvd = True

# time parameter setting
tidal_amplitude = 1.5  # galveston bay area normal tides.
tidal_period = 12.42*60*60
startime = 0.0
endtime = tidal_period*5
timestep = tidal_period/100
theta = 0.6

# equation parameter setting
gravity = 9.81
rhoair  = 1.225
rhowater = 1027
windx = 70.0
windy = -70.0


# bathymetry setting
#bathymetry = {"flat": 20.}
#bathymetry = {"vary": "2*x+3*y"}
#bathymetry = {"class": ['type1',400,600,'5.0','-3*((x-500)/100)**4 + 6*((x-500)/100)**2 + 2']}
bathymetry = {"class": ['type2', 2150, '-14.0/2150.0*x + 19.0', '5.0']}


# initial condition
initial_u = {"flat": (0., 0.)}
initial_eta = {"flat": 0.0}
#initial_u ={"vary": ('2*x+3*y', '3*y')}
#initial_eta = {"vary": "q0 * q1 * sp.cos(pi/100.0*x) "}

# boundary condition
# first way:
bcfile = "inlet_coarse_50_facet_region.xml"
boundary_u = {1: "freeslipyy", 2: "freeslipxx"} # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
boundary_eta = {3: "inlet special"} # "float/str"
# second way:
#bc_file = None
#boundary_u = {"Left": "freeslipyy", "Right": "freeslipyy", "Up": "freeslipxx", "Down": "freeslipxx"}# include sin cos ==> should add sympy symbol!!! # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
#boundary_eta = {}# "float/str"
