

from numpy import ndarray
from chaospy import E


def make_sto_ijk(pc_basis_str):

    n_modes = pc_basis_str.get("n_modes", False)
    orthogonal_pol = pc_basis_str.get("basis", False)
    joint_cdf = pc_basis_str.get("joint_cdf", False)

    if all([n_modes, orthogonal_pol, joint_cdf]) is False:
        raise TypeError("obtain NULL value from pc basis structure! Exit program...")

    sto_ijk_list = ndarray(shape=(n_modes, n_modes, n_modes), dtype=float)
    for i in range(n_modes):
        for j in range(n_modes):
            for k in range(n_modes):
                sto_ijk_list[i, j, k] = E(orthogonal_pol[i]*orthogonal_pol[j]*orthogonal_pol[k], joint_cdf)
    return sto_ijk_list
