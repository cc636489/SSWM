

from fenics import interpolate, Expression
from make_sto_modes import make_sto_modes
from random import random, seed


def make_initial_object_list(init_u_dict, init_eta_dict, pc_basis_str, u_function_space, eta_function_space,
                             randomness):

    n_modes = pc_basis_str.get("n_modes")

    # check if the keys is entered right.
    key_u = list(init_u_dict.keys())
    key_eta = list(init_eta_dict.keys())
    if len(key_u) != 1:
        raise TypeError("enter more than one initial condition for u! not allowed! Exit program...")
    if len(key_eta) != 1:
        raise TypeError("enter more than one initial condition for eta! not allowed! Exit program...")
    if key_u[0] != "flat" and key_u[0] != "vary":
        raise TypeError("don't have this kind of initial condition for u implemented! Exit program...")
    if key_eta[0] != "flat" and key_eta[0] != "vary":
        raise TypeError("don't have this kind of initial condition for eta implemented! Exit program...")

    # check if the values is entered right
    value_u = init_u_dict.get(key_u[0])
    if len(value_u) != 2:
        raise TypeError("enter the wrong length of initial condition for u! Tuple size should be 2! Exit program...")
    if not isinstance(value_u, tuple):
        raise TypeError("enter the wrong type initial conditions for u! Should be a tuple pair! Exit program...")
    if key_u[0] == "flat":
        if any([isinstance(value_u[ii], float) is False for ii in range(len(value_u))]):
            raise TypeError("enter the wrong type of initial condition for u_vel and v_vel! Should be float for both!")
    elif key_u[0] == "vary":
        if any([isinstance(value_u[ii], str) is False for ii in range(len(value_u))]):
            raise TypeError("enter the wrong type initial condition for u! Should be string for both! Exit program...")
    else:
        pass
    value_eta = init_eta_dict.get(key_eta[0])
    if key_eta[0] == "flat":
        if not isinstance(value_eta, float):
            raise TypeError("enter the wrong type initial condition for eta! Should be a float! Exit program...")
    elif key_eta[0] == "vary":
        if not isinstance(value_eta, str):
            raise TypeError("enter the wrong type initial condition for eta! Should be a string! Exit program...")
    else:
        pass

    # make the stochastic initial condition for u and eta.
    if key_u[0] == "flat":
        temp_u, temp_v = make_sto_modes(pc_basis_str, str(value_u[0]), str(value_u[1]))
    elif key_u[0] == "vary":
        temp_u, temp_v = make_sto_modes(pc_basis_str, value_u[0], value_u[1])
    else:
        raise TypeError("temp_u and temp_v has wrong length! Should be equal to n_modes! not allowed!")
    if len(temp_u) != n_modes or len(temp_v) != n_modes:
        raise TypeError("temp_u and temp_v has wrong length! Should be equal to n_modes! not allowed!")

    if key_eta[0] == "flat":
        temp_eta = make_sto_modes(pc_basis_str, str(value_eta))
    elif key_eta[0] == "vary":
        temp_eta = make_sto_modes(pc_basis_str, value_eta)
    else:
        raise TypeError("temp_eta has wrong length! Should be equal to n_modes! not allowed!")
    if len(temp_eta) != n_modes:
        raise TypeError("temp_eta has wrong length! Should be equal to n_modes! not allowed!")

    u_ic_list = []
    eta_ic_list = []
    for ii in range(n_modes):
        u_ic_list.append(temp_u[ii])
        u_ic_list.append(temp_v[ii])
        eta_ic_list.append(temp_eta[ii])

    if randomness:
        random_u = RandomUInitialConditions(element=u_function_space.ufl_element(), n_modes=n_modes)
        u_ic_expression = random_u
    else:
        u_ic_expression = Expression(u_ic_list, element=u_function_space.ufl_element())

    u_ic_function = interpolate(u_ic_expression, u_function_space)

    if n_modes == 1:
        eta_ic_function = interpolate(Expression(eta_ic_list[0], element=eta_function_space.ufl_element()),
                                      eta_function_space)
    else:
        eta_ic_function = interpolate(Expression(eta_ic_list, element=eta_function_space.ufl_element()),
                                      eta_function_space)

    return u_ic_function, eta_ic_function


class RandomUInitialConditions(Expression):

    def __init__(self, **kwargs):
        self.element = kwargs["element"]
        self.n_modes = kwargs["n_modes"]
        seed(self.n_modes)

    def eval(self, values, x):
        for i in range(self.n_modes * 2):
            values[i] = - 1e-6 + (random() - 0.1) / 0.9 * 2 * 1e-6

    def value_shape(self):
        return (self.n_modes * 2,)
