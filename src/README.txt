Plan of coding:

1. extend MakeDomain.py, enable Galveston Bay mesh input into fenics recognized version.

2. extend MakeBathy.py, enable Galveston Bay bathymetry input into fenics recognized version, put into ZERO mode, others are ZERO.

3*. extend Makeic.py, enable discretized initial condition input as ZERO mode, put every other mode ZERO.
    ==> Good to set it to zero everywhere.

4*. extend Makebc.py, enable discretized boundary condition input as ZERO mode, put every other mode ZERO. 
    ==> extend freeslip boundary condition(u.n=0) for curly boundaries. 
    ==> if none of these helps, then keep it "noslip".

5. build up MakeWind.py, calculated windstress at each node, at every timestep.
    a. reading in holland wind, x/y version, lat/lon version.
    b. calculate wind velocity at each node.
    c. generalized to every time step.
    ==> obtain |WindSpeed|*Wx and |WindSpeed|*Wy.
