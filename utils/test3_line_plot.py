from fenics import *
import numpy as np
import matplotlib.pyplot as plt

inputdir = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/'
testcase = 'test3_res_wiLES/'
commonfolder = '20181001_generalized/'
modestatus = 'ADCIRC_compare_with_same_mesh/'
# file = ['eta_500_100_300sec_adcirc0.csv', 'eta_500_100_300sec_sswm0.csv']
file = ['u_500_100_300sec_adcirc0.csv', 'u_500_100_300sec_sswm0.csv']
legend_name = ['eta(m)', 'u(m/s)', 'v(m/s)', 'umag(m/s)']
legend_index = 3

plt.figure()

for i in range(len(file)):

    ffield = inputdir + testcase + commonfolder + modestatus + file[i]
    time, field = [], []
    temp1, temp2 = 0, 0

    with open(ffield, 'r') as f:
        for j, line in enumerate(f):
            if j != 0:
                line = line.split(',')
                time.append(float(line[0]))
                if legend_index == 3:
                    temp1 = float(line[1])
                    temp2 = float(line[2])
                    field.append(np.sqrt(temp1**2+temp2**2))
                elif legend_index == 0 or legend_index == 1:
                    field.append(float(line[1]))
                elif legend_index == 2:
                    field.append(float(line[2]))

    print time
    print field

    plt.plot(time, field, '*-')

plt.legend(['ADCIRC', '2D-DSWM'], loc='upper center')
plt.xlim([min(time), max(time)])
plt.xticks(np.arange(min(time), (max(time)+0.001), step=(max(time)-min(time))/10))
plt.xlabel('time(sec)')
plt.ylim([-0.3, 0.7])
plt.ylabel(legend_name[legend_index])
plt.grid()
# plt.show()
plt.savefig(inputdir + testcase + commonfolder + modestatus + 'umag_comparison_for2_adcirc_sswm_500_100.png')

