
import numpy
from chaospy import Uniform, orth_ttr, J


def make_sto_basis(dist_name, sto_poly_deg, sto_poly_dim, co_eff):

    if dist_name == "uniform":
        if len(co_eff) / 2 != sto_poly_dim:
            raise TypeError("size of co_eff is not matching the dimension of polynomial chaos basis!")
        if any([co_eff[2 * z] > co_eff[2 * z + 1] for z in range(len(co_eff) // 2)]):
            raise TypeError("co_eff value should be in order as lower-high-lower-high, invalid number entering!")
    elif dist_name == "gaussian":
        if len(co_eff)/2 != sto_poly_dim:
            raise TypeError("size of co_eff is not matching the dimension of polynomial chaos basis!")
    else:
        raise TypeError("Only implement uniform and gaussian distribution! dist_name not found!")

    if dist_name == "uniform":
        ker_cis = [Uniform(co_eff[2 * z], co_eff[2 * z + 1]) for z in range(len(co_eff) // 2)]
        cp_string = "J("
        for ii in range(len(ker_cis)):
            cp_string = cp_string + "ker_cis[" + str(ii)+"]"
            if ii < len(ker_cis)-1:
                cp_string = cp_string + ","
        cp_string = cp_string + ")"
        joint_cdf = eval(cp_string)
        joint_pdf = 1.0
        for ii in range(len(ker_cis)):
            joint_pdf = joint_pdf / (co_eff[2 * ii + 1] - co_eff[2 * ii])
        orthogonal_pol = orth_ttr(sto_poly_deg, joint_cdf, normed=True)
        n_modes = len(orthogonal_pol)
    elif dist_name == "gaussian":
        raise TypeError("Not implement Gaussian distribution yet! Exit program...")
    else:
        raise TypeError("Only implement uniform and gaussian distribution! dist_name not found!")
    pc_basis_str = {"name": dist_name, "dim": sto_poly_dim, "order": sto_poly_deg, "n_modes": n_modes, "co_eff": co_eff,
                    "basis": orthogonal_pol, "joint_cdf": joint_cdf, "joint_pdf": joint_pdf}

    return pc_basis_str
