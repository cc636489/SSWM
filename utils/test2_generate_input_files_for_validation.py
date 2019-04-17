
import numpy as np
import subprocess


output_dir1 = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"
output_dir2 = "test2_truth/"
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
input_file = "input_generalized_testcase_stochastic_2_done.py"

for q0 in np.linspace(-1.2, 1.2, num = 20):
    for q1 in np.linspace(-2.0, 2.0, num = 20):
        output_file = input_file[:27] + "deterministic_2_" + str(round(q0, 2)) + "_" + str(round(q1, 2)) + ".py"
        subprocess.call(["cp " + input_dir + input_file + " " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's#output_dir = \"/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/\"#output_dir = \"/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/test2_truth/\"#' " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's/output_str = \"test2_stochastic_\"/output_str = \"test2_deterministic_" + str(round(q0, 2)) + "_" + str(round(q1, 2)) + "_\"/' " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's/sto_poly_deg = 4/sto_poly_deg = 0/' " + output_dir1 + output_dir2 + output_file], shell = True)
        string = str(0.1 * q0 * q1)
        print "initial_eta = " + string + "*cos(pi/100.0*x)"
        subprocess.call(["sed -i 's#initial_eta = {\"vary\": \"0.1\\*q0\\*q1\\*sp\\.cos(pi/100.0\\*x)\"}#initial_eta = {\"vary\": \"" + string + "\\*sp\\.cos(pi/100.0\\*x)\"}#' " + output_dir1 + output_dir2 + output_file], shell = True)

