
import numpy as np
import subprocess


output_dir1 = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"
output_dir2 = "test3_truth/"
#input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
input_dir = "./"
input_file = "input_generalized_testcase_stochastic_3_done_for_input_generation.py"

for q0 in np.linspace(0.8, 1.2, num = 20):
    for q1 in np.linspace(0.9, 1.1, num = 20):
        output_file = input_file[:27] + "deterministic_3_" + str(round(q0, 2)) + "_" + str(round(q1, 2)) + ".py"
        subprocess.call(["cp " + input_dir + input_file + " " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's#output_dir = \"/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/\"#output_dir = \"/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/test3_truth/\"#' " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's/output_str = \"test3_stochastic_\"/output_str = \"test3_deterministic_" + str(round(q0, 2)) + "_" + str(round(q1, 2)) + "_\"/' " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's/sto_poly_deg = 3/sto_poly_deg = 0/' " + output_dir1 + output_dir2 + output_file], shell = True)
        string1 = str(-3.0 * q0 + 6.0 * q1 + 2.0)
        string2 = str(-3.0 * q0)
        string3 = str(6.0 * q1)
        print "bathy = " + string1 + ', ' + string2 + ' * ((x-500.0)/100.0)**4 + ' + string3 + ' * ((x-500.0)/100.0)**2 + 2.0'
        subprocess.call(["sed -i 's/string1 = \"-3.0 \\* q0 + 6.0 \\* q1 + 2.0\"/string1 = \"" + string1 + "\"/g' " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's/string2 = \"-3.0 \\* q0\"/string2 = \"" + string2 + "\"/g' " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's/string3 = \"6.0 \\* q1\"/string3 = \"" + string3 + "\"/g' " + output_dir1 + output_dir2 + output_file], shell = True)

