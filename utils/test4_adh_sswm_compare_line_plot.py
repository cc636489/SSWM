from fenics import *
import numpy as np
import matplotlib.pyplot as plt

inputdir = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/'
testcase = 'test8_SUPG_SamgLily_crosswind/'
commonfolder = 'results/'
# modestatus = 'adcirc_result_comparison/'
file = ['adh-all-250-0.csv', 'sswm-eta-250-0.csv']
# file = ['uv_adh_entrance-250.0.csv', 'uv_sswm_entrance-250.0.csv']
legend_name = ['eta(m)', 'u(m/s)']
legend_index = 0

plt.figure(figsize=[8, 6])
for i in range(len(file)):

    ffield = inputdir + testcase + commonfolder + file[i]
    time, field = [], []

    with open(ffield, 'r') as f:
        for j, line in enumerate(f):
            if j != 0:
                line = line.split(',')
                time.append(float(line[0]))
                field.append(float(line[1]))

    print time
    print field

    time = [i/3600/24 for i in time]
    plt.plot(time, field, '*-')

plt.legend(['AdH', '2D_DSWM'], loc='upper center')
plt.xlim([min(time), max(time)])
plt.xticks(np.arange(min(time), (max(time)+0.001), step=(max(time)-min(time))/10))
plt.xlabel('time(days)')
plt.ylim([-1.0, 1.0])
plt.ylabel(legend_name[legend_index])
plt.grid()
# plt.show()
plt.savefig(inputdir + testcase + commonfolder + 'eta_comparison_for2_adh_sswm_250_0.png')
