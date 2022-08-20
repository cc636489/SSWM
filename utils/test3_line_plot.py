from fenics import *
import matplotlib
matplotlib.use('Agg')
font = {'family': 'normal', 'size': 22}
matplotlib.rc('font', **font)
import numpy as np
import matplotlib.pyplot as plt

inputdir = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/'
testcase = 'test3_res/'
modestatus = 'SSWM_results/'

point_x = 250 
point_y = 100
ylabel_name = ['eta(m)', 'u(m/s)', 'v(m/s)', 'umag(m/s)']
save_name = ['eta', 'u', 'v', 'umag']
index = 2 
file_name = ['eta', 'u']
file_index = 1
num = 6

file = [file_name[file_index]+'_'+str(point_x)+'_'+str(point_y)+'_300sec_sswm.csv', file_name[file_index]+'_'+str(point_x)+'_'+str(point_y)+'_300sec_adcirc.csv']

plt.subplots(figsize=(11.5, 7))

for i in range(len(file)):

    ffield = inputdir + testcase + modestatus + file[i]
    time, field = [], []
    temp1, temp2 = 0, 0

    with open(ffield, 'r') as f:
        for j, line in enumerate(f):
            if j != 0:
                line = line.split(',')
                if index == 0:
                    field.append(float(line[0]))
                    time.append(float(line[2]))
                elif index == 1:
                    field.append(float(line[0]))
                    time.append(float(line[4]))
                elif index == 2:
                    field.append(float(line[1]))
                    time.append(float(line[4]))
                elif index == 3:
                    temp1 = float(line[0])
                    temp2 = float(line[1])
                    field.append(np.sqrt(temp1**2+temp2**2))
                    time.append(float(line[4]))
    if i == 0:
        plt.plot(time, field, '-', color='mediumseagreen', linewidth=4)
        plt.xlim([min(time), max(time)])
        plt.xticks(np.linspace(min(time), max(time), num))
    else:
        plt.plot(time[:300], field[:300], '*', color='darkorange', markersize=12)
        plt.ylim([min(field[:300])*1.5, max(field[:300])*1.5])
        
plt.legend(['DSWM', 'ADCIRC'], loc='upper left')
plt.xlabel('time(seconds)')
plt.ylabel(ylabel_name[index])
plt.title('Comparison at point (' + str(point_x) + ',' + str(point_y) + ')')
plt.grid()
plt.savefig(inputdir + testcase + modestatus + save_name[index] + '_comparison_for2_adcirc_sswm_' + str(point_x) + '_' + str(point_y) + '.pdf', bbox_inches = 'tight')

