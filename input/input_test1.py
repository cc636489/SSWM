##################################################################
# Input parameter
##################################################################
# output file name string
outputdir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/test1_res/"
outputstr = "fem_verified_test1_"

# stochastic input
nsample = 20
xytickstep = 4
xyticknumber = 6
scatterplotsize = 5
#only implement with only one testnode at one time.
testnodeX = [25.]
testnodeY = [25.]
#zmax = 1.3
#zmin = -8.3

# stochastic basis
distname = "uniform"  # "uniform"  or "gaussian"   ==> only enable "uniform" mode at present.
sto_poly_deg = 4 # polynomial chaos order is 2.
# the dimension and coefficient should be paired and dim!=0  and coefficient!=Null
sto_poly_dim = 2 # use "q0","q1","q2", ....
coefficient = [-1.2, 1.2, -2, 2] # lower1/upper1--lower2/upper2--...

# horizontal domain setting
# first way: simple domain built in
domain = {"rectangle": (0., 0., 100., 50., 20, 10)}
# second way: complex domain import: domain = {"importfile": "inlet.xml"}

# stochastic coefficient # if contains sin(), cos(), should use sympy sin and sympy cos!!!!!!!!!!
sto_viscosity = "0.0"
sto_friction = "0.0"
sto_force_x = "0.0"# can be time dependent
sto_force_y = "0.0"# can be time dependent

# terms control
include_viscosity = True
include_advection = True
linear_divergence = False

# time parameter setting
startime = 0.0
endtime = 25.0
timestep = 0.5
theta = 0.5

# equation parameter setting
gravity = 9.81

# bathymetry setting
bathymetry = {"flat": 20.}
#bathymetry = {"vary": "2*x+3*y"}
#bathymetry = {"class": ['type1',400,600,'5.0','-3*((x-500)/100)**4 + 6*((x-500)/100)**2 + 2']}

# initial condition
initial_u = {"flat": (0., 0.)}
#initial_eta = {"flat": 0.0}
#initial_u ={"vary": ('2*x+3*y', '3*y')}
initial_eta = {"vary": "q0 * q1 * sp.cos(pi/100.0*x) "}

# boundary condition
# first way:
#bc_file = "inlet_facet_region.xml"
#boundary_u = {1: ("sp.sin(pi*t/5)*x*(50.0-y)/625.0", "0.0"), 2: "freeslipxx"} # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
#boundary_eta = {3: 0.0} # "float/str"
# second way:
bcfile = None
boundary_u = {"Left": "freeslipyy", "Right": "freeslipyy", "Up": "freeslipxx", "Down": "freeslipxx"}# include sin cos ==> should add sympy symbol!!! # "freeslipxx/freeslipyy/noslip/(str, tuple)/(float, tuple)"
boundary_eta = {}# "float/str"
