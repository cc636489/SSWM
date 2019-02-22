2D-SSWM: 2D Stochastic Shallow Water Model in Python
================================================================

2D-SSWM is a 2D shallow water model, used to study the uncertainty of surface elevation due to the uncertain model input.
The model has adopted stochastic finite element method (SFEM), incorporating Streamlined Upwind Petrov Galerkin (SUPG) term
and crosswind diffusion term for stabilization purpose.

**SFEM**, in general, is a combination of the following two methods:
-   spectrum method in random space
-   Continuous Galerkin(CG) finite element method in physical space.

**SUPG** is added to the model in order to stabilize convection-dominated flow in shallow water problems.
-   It's equivalent to add diffusion in the streamlined direction.

**Crosswind diffusion** term is added to the model in order to stabilize spurious oscillation due to large gradients and sharp fronts.
-   It's equivalent to add diffusion in the large gradient direction.

**When** to use it:
-   interested in the water surface elevation under the wind forces, like hurricanes.
-   want to study the combined effect of the uncertain model inputs on the model output. The uncertain model inputs include bathymetry, initial conditions,
time-varying boundary conditions, wind drag coefficient and bottom friction coefficient.
-   want to know the mean field and the variance map of surface elevation(or water velocity) due to uncertain model inputs.
-   want to know the probability density function of surface elevation(or water velocity) at one node due to uncertain model inputs.
-   want to know the auto-correlation(or cross-correlation) of the time series of surface elevation at any two nodes.
-   want to obtain a highly-accurate surrogate surface at any nodes.

Dependencies:
-------------

-   Fenics(2017.1.0)
-   chaospy(2.2.3)
-   numpy(1.13.3)
-   scipy(0.19.1)
-   sympy(1.1.1)
-   netCDF4(1.4.0)

Inputs:
------------

model input template:
-   ./input/input_generalized.py

model input files:
-   mesh file (in .xml format)
-   bathymetry file (in netCDF4 format)
-   boundary file (in .xml format)
-   Hurricane best tack input (in .txt format)

model input options:
-   turn on/off DEBUG mode.
-   turn on/off SUPG terms.
-   turn on/off crosswind terms.
-   turn on/off large eddy simulation Smagorinsky model.

Outputs:
--------

model output files:
-   surface elevations in each stochastic modes (in .h5 and/or .vtu format)
-   water velocity in each stochastic modes (in .h5 and/or .vtu format)


Execution:
----------

-   download the repository.
-   run all the test cases `python -v ./tests/test_*.py`.
-   prepare model input files.
-   modify model parameters in ../input/input_generalized.py
-   run the program `python driver.py`.

Model Verification:
-------------------

Verification of the program is conducted by comparing to a well-known shallow water model ADCIRC. The following is the
comparison of surface elevation between ADCIRC and SSWM(i.e. **mean surface elevation** with stochastic order setting
to zero)
![Alter Text](https://github.com/cc636489/research/blob/master/doc/IKE_comparison_u_eta_with_atmos_pre_faster_animation.gif)
