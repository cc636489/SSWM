
import numpy as np
import subprocess


output_dir1 = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"
output_dir2 = "HARVEY_truth/"
input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
input_file = "input_generalized_testcase_stochastic_Gulf_wind_HARVEY.py"

for q0 in np.linspace(0.8, 1.2, num = 50):
    #for q1 in np.linspace(0.9, 1.1, num = 10):
        #output_file = input_file[:27] + "deterministic_4_" + str(round(q0, 2)) + "_" + str(round(q1, 2)) + ".py"
        output_file = input_file[:27] + "deterministic_HARVEY_" + str(round(q0, 4)) + ".py"
        subprocess.call(["cp " + input_dir + input_file + " " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's#output_dir = \"/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/\"#output_dir = \"/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/HARVEY_truth/\"#' " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's/output_str = \"gulf_winds_harvey_stochastic_4.0_2.0_\"/output_str = \"gulf_winds_harvey_deterministic_" + str(round(q0, 4)) + "_\"/' " + output_dir1 + output_dir2 + output_file], shell = True)
        subprocess.call(["sed -i 's/sto_poly_deg = 1/sto_poly_deg = 0/' " + output_dir1 + output_dir2 + output_file], shell = True)
        #string = str(0.001 * q0 * q1)
        #print "sto_windDrag = " + string
        string = str(q0)
        print "sto_windDrag = " + string
        subprocess.call(["sed -i 's/sto_windDrag = \"q0\"/sto_windDrag = \"" + string + "\"/' " + output_dir1 + output_dir2 + output_file], shell = True)
        #subprocess.call(["sed -i 's/tidal_amplitude = \"0.2\\*q0\"/tidal_amplitude = " + string + "/g' " + output_dir1 + output_dir2 + output_file], shell = True)
        #subprocess.call(["sed -i 's/M2 special stochastic/M2 special/g' " + output_dir1 + output_dir2 + output_file], shell = True)

