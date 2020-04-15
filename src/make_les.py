

from fenics import Function, grad, TestFunction, TrialFunction, sqrt, \
    inner, LinearVariationalProblem, LinearVariationalSolver, CellVolume, dx


class LES(object):

    def __init__(self, v, u, s_coefficient):
        self.V = v
        self.u = u
        self.cs = s_coefficient
        self.eddy_viscosity = Function(v)

    def _strain_rate_tensor(self):
        ss = 0.5 * (grad(self.u) + grad(self.u).T)
        return ss

    def _eddy_viscosity_eqn(self):

        dim = len(self.u)
        ww = TestFunction(self.V)
        qq = TrialFunction(self.V)

        cell_vol = CellVolume(self.V.mesh())
        filter_width = cell_vol ** (1.0 / dim)

        ss = self._strain_rate_tensor()
        second_invariant = 0.0
        for ii in range(0, dim):
            for jj in range(0, dim):
                second_invariant += 2.0 * (ss[ii, jj] ** 2)

        second_invariant = sqrt(second_invariant)
        forms = (self.cs * filter_width) ** 2 * second_invariant

        l_h_s = inner(ww, qq) * dx
        r_h_s = inner(ww, forms) * dx

        return l_h_s, r_h_s

    def solve(self):

        les_lhs, les_rhs = self._eddy_viscosity_eqn()
        eddy_viscosity_problem = LinearVariationalProblem(les_lhs, les_rhs, self.eddy_viscosity, bcs=[])
        _solver = LinearVariationalSolver(eddy_viscosity_problem)
        _solver.parameters["linear_solver"] = "lu"
        _solver.parameters["symmetric"] = True
        _solver.parameters["lu_solver"]["reuse_factorization"] = True
        _solver.solve()

        return self.eddy_viscosity
