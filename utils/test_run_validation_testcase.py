
import numpy as np
import subprocess


input_dir1 = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"
input_dir2 = "HARVEY_truth/"
input_file = "input_generalized_testcase_"
target_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
exec_file = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/src/driver.py"

for q0 in np.linspace(0.8244898, 0.84897959, num = 4):
    #for q1 in np.linspace(1.0, 2.0, num = 20):
        #output_input_file = input_file + "deterministic_4_" + str(round(q0, 2)) + "_" + str(round(q1, 2)) + ".py"
        output_input_file = input_file + "deterministic_HARVEY_" + str(round(q0, 4)) + ".py"
        subprocess.call(["cp " + input_dir1 + input_dir2 + output_input_file + " " + target_dir + "input_generalized.py"], shell = True)
        subprocess.call(["python " + exec_file], shell = True)
