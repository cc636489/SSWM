
import numpy as np
import subprocess


input_dir1 = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test_stochastic/"
input_dir2 = "test3_truth/"
input_file = "input_generalized_testcase_stochastic_3_done.py"
target_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/input/"
exec_file = "/workspace/Documentation/Research_Doc/SFEM_Doc/7-SSWM-github/src/driver.py"

for q0 in np.linspace(0.8, 1.2, num = 20):
    for q1 in np.linspace(0.9, 1.1, num = 20):
        output_input_file = input_file[:27] + "deterministic_3_" + str(round(q0, 2)) + "_" + str(round(q1, 2)) + ".py"
        subprocess.call(["cp " + input_dir1 + input_dir2 + output_input_file + " " + target_dir + "input_generalized.py"], shell = True)
        subprocess.call(["python " + exec_file], shell = True)
