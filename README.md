2D-SSWM: 2D Stochastic Shallow Water Model in Python
================================================================

2D-SSWM is a 2D stochastic shallow water model, used to study the uncertainty of surface elevation/water velocity due to the uncertain model input.
The model has adopted stochastic finite element method (SFEM), introduced Stochatic Streamlined Upwind Petrov Galerkin (SSUPG) 
and Stochastic Crosswind Diffusion(SCD) for stabilization purpose.

**SFEM**, in general, is a combination of the following two methods:
-   Stochastic finite element method in probability space
-   Continuous Galerkin(CG) finite element method in physical space.

**SSUPG** is introduced to the model in order to stabilize convection-dominated flow in each modes in stochastic shallow water problems.
-   It's equivalent to add diffusion in the streamlined direction.

**SCD** is introduced to the model in order to stabilize spurious oscillation due to large gradients in each modes and sharp fronts.
-   It's equivalent to add diffusion in the large gradient direction.

**When** to use it:
-   interested in the water surface elevation under the wind forces, like hurricanes.
-   want to study the combined effect of the uncertain model inputs on the model output. The uncertain model inputs include bathymetry, initial conditions,
time-varying boundary conditions, wind drag coefficient and bottom friction coefficient.
-   want to know the mean field and the variance map of surface elevation(or water velocity) due to uncertain model inputs.
-   want to know the probability density function of surface elevation(or water velocity) at one location due to uncertain model inputs.
-   want to know the auto-correlation(or cross-correlation) of the time series of surface elevation at any two nodes.
-   want to obtain a highly-accurate surrogate surface at any locations.

Dependencies:
-------------

-   Fenics(2018.1.0)
-   chaospy(2.2.4)
-   numpy(1.15.4)
-   scipy(1.2.1)
-   sympy(1.10.1)
-   netCDF4(1.5.1.2)
-   python(3.7.12)

The best way to get all the dependencies ready on your local machine is through conda. Then use for example `conda activate fenics2018project` to activate your conda environment.

Inputs:
------------

model input template:
-   ./input/input_generalized.py

model input files:
-   mesh file (in .xml format)
-   bathymetry file (in netCDF4 format)
-   boundary file (in .xml format)
-   Hurricane best track input (in .txt format)

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

-   clone the repository.
-   prepare model input files.
-   modify model parameters in ../input/input_generalized.py
-   run the program `python driver.py`.

SSWM Model Validation( Deterministic Part )
---------------------------------------------

Validation of the deterministic part of SSWM is conducted by comparing to:
-   analytical solution with convergence analysis.
-   well-established numerical model
-   experimental data

The analytical solution with convergence analysis is opted out here. please find the reference [1] for more information.

The following is the
comparison of surface elevation between ADCIRC and SSWM(i.e. **mean surface elevation** with stochastic order setting
to zero)
<p align="center">
<img src="https://github.com/cc636489/research/blob/master/doc/IKE_comparison_u_eta_with_atmos_pre_faster_animation.gif">
</p>

The following is the velocity profile comparison against experimental data.
<p align="center">
<img src="https://github.com/cc636489/research/blob/master/doc/Experimental_comparison.png">
</p>


SSWM Model Validation( Stochastic Part )
---------------------------------------------

Validation of the stochastic part of SSWM is conducted by comparing pdfs against its monto-carlo counterparts. The uncertain resources we considered here includes:
-   uncertain initial condition (slosh test).
<p align="center">
<img src="https://github.com/cc636489/research/blob/master/doc/PDF_comparison_slosh.png">
</p>

-   uncertain bathymetry condition (hump test).
<p align="center">
<img src="https://github.com/cc636489/research/blob/master/doc/PDF_comparison_hump.png">
</p>

-   uncertain boundary condition (inlet test).
<p align="center">
<img src="https://github.com/cc636489/research/blob/master/doc/PDF_comparison_Inlet.png">
</p>

-   uncertain model parameter (hurricane tests).
<p align="center">
<img src="https://github.com/cc636489/research/blob/master/doc/PDF_comparison_hurricane_Harvey.png">
</p>


Reference
----------

[1] Chen C, Dawson C, Valseth E. Cross-mode stabilized stochastic shallow water systems using stochastic finite element methods[J]. Computer Methods in Applied Mechanics and Engineering, 2023, 405: 115873. (https://doi.org/10.1016/j.cma.2022.115873)
