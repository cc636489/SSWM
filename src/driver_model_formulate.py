

from model_formulate import DetModelFormulate, StoModelFormulate


def model_formulate_driver(inputs, initiate):

    if initiate.basis_str.get('n_modes') == 1:

        formulate = DetModelFormulate(inputs, initiate)

        if inputs.include_convection:
            formulate.add_convection()

        if inputs.include_viscosity:
            formulate.add_viscosity()

        if inputs.include_bottom_stress:
            formulate.add_bottom_stress()

        if inputs.include_wind_stress:
            formulate.add_wind_stress()

        if inputs.include_atmospheric_pressure:
            formulate.add_atmospheric_pressure()

        if inputs.include_supg:
            formulate.add_su_pg()

        if inputs.include_crosswind:
            formulate.add_crosswind()

    else:

        formulate = StoModelFormulate(inputs, initiate)

        if inputs.include_convection:
            formulate.add_convection()

        if inputs.include_viscosity:
            formulate.add_viscosity()

        if inputs.include_bottom_stress:
            formulate.add_bottom_stress()

        if inputs.include_wind_stress:
            formulate.add_wind_stress()

        if inputs.include_supg:
            formulate.add_su_pg()

        if inputs.include_crosswind:
            formulate.add_crosswind()

        if inputs.include_atmospheric_pressure:
            formulate.add_atmospheric_pressure()

    return formulate
