
import numpy as np
import subprocess


output_dir1 = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"
output_dir2 = "test4_truth/"
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
input_file = "input_generalized_testcase_stochastic_4_done.py"

for q0 in np.linspace(1.0, 2.0, num = 50):
    #for q1 in np.linspace(0.9, 1.1, num = 10):
        #output_file = input_file[:27] + "deterministic_4_" + str(round(q0, 2)) + "_" + str(round(q1, 2)) + ".py"
        output_file = input_file[:27] + "deterministic_4_" + str(round(q0, 2)) + ".py"
        subprocess.call(["cp " + input_dir + input_file + " " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's#output_dir = \"/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/\"#output_dir = \"/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/test4_truth/\"#' " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's/output_str = \"test4_stochastic_\"/output_str = \"test4_deterministic_" + str(round(q0, 2)) + "_\"/' " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's/sto_poly_deg = 3/sto_poly_deg = 0/' " + output_dir1 + output_dir2 + output_file], shell = True)
        #string = str(0.001 * q0 * q1)
        #print "sto_windDrag = " + string
        string = str(0.2 * q0)
        print "tidal_amplitude = " + string
        #subprocess.call(["sed -i 's/sto_windDrag = \"0.001 \\* q0 \\* q1\"/sto_windDrag = \"" + string + "\"/' " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's/tidal_amplitude = \"0.2\\*q0\"/tidal_amplitude = " + string + "/g' " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's/M2 special stochastic/M2 special/g' " + output_dir1 + output_dir2 + output_file], shell = True)

