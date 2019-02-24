

from fenics import MeshFunction, FacetFunction, SubDomain, near, DirichletBC, Expression, Constant, pi
from make_sto_modes import make_sto_modes
from numpy import zeros, pi


def make_boundary_object_list(boundary_u_dict, boundary_eta_dict, pc_basis_str, bc_file_name, mesh, domain,
                              u_function_space, eta_function_space, current_t, tidal_amplitude, tidal_period):

    # check input boundary condition keys.
    n_modes = pc_basis_str.get("n_modes")
    n_u_keys = len(boundary_u_dict.keys())
    n_eta_keys = len(boundary_eta_dict.keys())
    if not (isinstance(bc_file_name, str) or bc_file_name is None):
        raise TypeError("wrong type bc_file_name! Exit Program...")
    if all([isinstance(boundary_u_dict.keys()[i], int) for i in range(n_u_keys)]) is False:
        raise TypeError("wrong key type of boundary u! Exit Program...")
    if all([isinstance(boundary_eta_dict.keys()[i], int) for i in range(n_eta_keys)]) is False:
        raise TypeError("wrong key type of boundary eta! Exit Program...")

    # initialize returned bc object list.
    u_bc_list = []
    eta_bc_list = []
    u_list_expression = []
    eta_list_expression = []
    u_time_dependent = False
    eta_time_dependent = False

    # make boundaries sub_domain and mark it sequentially.
    if isinstance(bc_file_name, str):
        boundaries = MeshFunction("size_t", mesh, bc_file_name)
    elif bc_file_name is None:
        boundaries = FacetFunction("size_t", mesh)

        # this 4 class only used for strong imposed non-penetration bc.==> maybe not
        class Left(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[0], domain.get("rectangle")[0])

        class Right(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[0], domain.get("rectangle")[2])

        class Up(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[1], domain.get("rectangle")[3])

        class Down(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[1], domain.get("rectangle")[1])

        # fenics build-in mesh don't have boundary marked, marked here sequentially, 1, 2, 3, 4...
        boundaries.set_all(0)
        for i, item in enumerate([Left, Right, Up, Down]):
            temp_sub_domain = item()
            temp_sub_domain.mark(boundaries, i+1)
    else:
        raise TypeError("bc_file_name is not allowed. Check input variable.")

    # make boundary object list for u.
    for key in boundary_u_dict.keys():
        value = boundary_u_dict.get(key)
        if value == "no_slip":
            u_list = list(zeros(n_modes*2))
            u_bc_list.append(DirichletBC(u_function_space, u_list, boundaries, key))
        elif value == "free_slip_in_x":
            if n_modes == 1:
                u_bc_list.append(DirichletBC(u_function_space.sub(1), Constant(0), boundaries, key))
            else:
                for j in range(n_modes):
                    u_bc_list.append(DirichletBC(u_function_space.sub(j).sub(1), Constant(0), boundaries, key))
        elif value == "free_slip_in_y":
            if n_modes == 1:
                u_bc_list.append(DirichletBC(u_function_space.sub(0), Constant(0), boundaries, key))
            else:
                for j in range(n_modes):
                    u_bc_list.append(DirichletBC(u_function_space.sub(j).sub(0), Constant(0), boundaries, key))
        elif isinstance(value, tuple) and len(value) == 2 and all([isinstance(value[j], float) for j in range(2)]):
            u_list = list(value) + list(zeros((n_modes-1)*2))
            u_bc_list.append(DirichletBC(u_function_space, u_list, boundaries, key))
        elif isinstance(value, tuple) and len(value) == 2 and all([isinstance(value[j], str) for j in range(2)]):
            temp_u_list, temp_v_list = make_sto_modes(pc_basis_str, value[0], value[1])
            u_list = []
            for j in range(n_modes):
                u_list.append(temp_u_list[j])
                u_list.append(temp_v_list[j])
            if any(["t" in str_list for str_list in u_list]):
                u_time_dependent = True
                u_list_expression = Expression(u_list, element=u_function_space.ufl_element(), t=current_t)
                u_bc_list.append(DirichletBC(u_function_space, u_list_expression, boundaries, key))
            else:
                u_list_expression = Expression(u_list, element=u_function_space.ufl_element())
                u_bc_list.append(DirichletBC(u_function_space, u_list_expression, boundaries, key))
        else:
            raise TypeError("wrong boundary u type.")

    # make boundary object list for eta.
    for key in boundary_eta_dict.keys():
        value = boundary_eta_dict.get(key)
        if isinstance(value, float):
            eta_list = [value] + list(zeros(n_modes-1))
            if n_modes == 1:
                eta_bc_list.append(DirichletBC(eta_function_space, eta_list[0], boundaries, key))
            else:
                eta_bc_list.append(DirichletBC(eta_function_space, eta_list, boundaries, key))
        elif isinstance(value, str) and value != "M2 special":
            eta_list = make_sto_modes(pc_basis_str, value)
            if any(["t" in str_list for str_list in eta_list]):
                eta_time_dependent = True
                if n_modes == 1:
                    eta_list_expression = Expression(eta_list[0], element=eta_function_space.ufl_element(), t=current_t)
                else:
                    eta_list_expression = Expression(eta_list, element=eta_function_space.ufl_element(), t=current_t)
            else:
                if n_modes == 1:
                    eta_list_expression = Expression(eta_list[0], element=eta_function_space.ufl_element())
                else:
                    eta_list_expression = Expression(eta_list, element=eta_function_space.ufl_element())
            eta_bc_list.append(DirichletBC(eta_function_space, eta_list_expression, boundaries, key))
        elif value == "M2 special":
            eta_time_dependent = True
            if n_modes == 1:
                eta_list_expression = Expression("amp*sin(omega*t)", element=eta_function_space.ufl_element(),
                                                 t=Constant(0), amp=tidal_amplitude, omega=2 * pi / tidal_period)
            else:
                temp_list = ["amp*sin(omega*t)"]
                for mode in range(n_modes-1):
                    temp_list.append("0.0")
                eta_list_expression = Expression(temp_list, element=eta_function_space.ufl_element(),
                                                 t=current_t, amp=tidal_amplitude, omega=2 * pi / tidal_period)
            # if inlet test case, boundary number should be 2;
            # if generalized gulf domain, the boundary number should be 3.
            eta_bc_list.append(DirichletBC(eta_function_space, eta_list_expression, boundaries, key))
        else:
            raise TypeError("enter wrong boundary eta type.")

    return u_bc_list, eta_bc_list, u_time_dependent, eta_time_dependent, u_list_expression, eta_list_expression

