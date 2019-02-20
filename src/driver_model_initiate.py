""" function to import mesh, set up bathymetry, initial condition and boundary condition """

from model_initiate import ModelInitiate


def model_initiate_driver(inputs):

    initiate = ModelInitiate(inputs)

    initiate.make_stochastic_basis()

    initiate.make_mesh()

    initiate.make_function_space()

    initiate.make_function()

    initiate.make_bath()

    initiate.make_ic()

    initiate.make_bc()

    if inputs.include_wind_stress and inputs.include_atmospheric_pressure:
        initiate.make_wind()

    initiate.make_stochastic_parameter()

    return initiate
