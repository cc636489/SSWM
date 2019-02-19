from makestobasis import MakeStoBasis
import chaospy as cp


def test_orthognality():

    distname = "uniform"
    sto_poly_deg = 2
    sto_poly_dim = 2
    coefficient = [0.18, 0.22, -0.001, 0.001]
    basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
    orth = basis_struct.get('basis')
    nmodes = basis_struct.get('nmodes')
    cdf = basis_struct.get('jointcdf')
    for i in range(nmodes):
        for j in range(nmodes):
            if i == j:
                if abs(cp.E(orth[i]*orth[j], cdf)-1.0) > 1.e-10:
                    print "Not satisfy orthogonality condition for basis " + str(i) + " and basis " + str(j)
                    return 0
            else:
                if abs(cp.E(orth[i]*orth[j], cdf)) > 1.e-10:
                    print "Not satisfy orthogonality condition for basis " + str(i) + " and basis " + str(j)
                    return 0

    return 1


def main():

    red = '\033[91m'
    green = '\033[92m'
    blue = '\033[94m'
    bold = '\033[1m'

    flag = test_orthognality()
    if flag:
        print blue+bold+"Test: "+bold+green+"PASS"
    else:
        print blue+bold+"Test: "+bold+red+"FAIL"


if __name__ == '__main__':
    main()