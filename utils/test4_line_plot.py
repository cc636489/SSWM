from fenics import *
import matplotlib
matplotlib.use('Agg')
font = {'family' : 'normal', 'size'   : 22}
matplotlib.rc('font', **font)
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

inputdir = '/Users/chenchen/gloria/4-NS-results-and-tests/' \
           'test4_res/adh_sswm_comparison_SUPG_SamgLily_crosswind/line_plot_results/'

point_x = 450
point_y = 0
ylabel_name = ['eta(m)', 'u(m/s)', 'v(m/s)', 'umag(m/s)']
save_name = ['eta', 'u', 'v', 'umag']
index = 1
file_name = ['eta', 'u']
file_index = 1
scale1 = 1.75 
scale2 = 1.75
num = 6

file = ['sswm-'+file_name[file_index]+'-'+str(point_x)+'-'+str(point_y)+'.csv', 'adh-'+file_name[file_index]+'-'+str(point_x)+'-'+str(point_y)+'.csv']

plt.subplots(figsize=(11.5, 7))

for i in range(len(file)):

    ffield = inputdir + file[i]
    time, field = [], []
    temp1, temp2 = 0, 0

    with open(ffield, 'r') as f:
        for j, line in enumerate(f):
            if j != 0:
                line = line.split(',')
                if index == 0:
                    field.append(float(line[1]))
                    time.append(float(line[0])/24./3600.)
                elif index == 1:
                    field.append(float(line[0]))
                    time.append(float(line[7])/24./3600.)
                elif index == 2:
                    field.append(float(line[2]))
                    time.append(float(line[0])/24./3600.)
                elif index == 3:
                    temp1 = float(line[1])
                    temp2 = float(line[2])
                    field.append(np.sqrt(temp1**2+temp2**2))
                    time.append(float(line[0])/24./3600.)

    if i == 0:
        plt.plot(time, field, '-', color='mediumseagreen', linewidth=4)
        plt.xlim([min(time), max(time)])
        plt.xticks(np.linspace(min(time), max(time), num))
        plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
    else:
        plt.plot(time[::4], field[::4], '*', color='darkorange', markersize=12)
        plt.ylim([min(field)*scale1, max(field)*scale2])
        plt.yticks(np.linspace(min(field)*scale1, max(field)*scale2, num))
        plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2f'))
        
plt.legend(['DSWM', 'ADH'], loc='upper right')
plt.xlabel('time(days)')
plt.ylabel(ylabel_name[index])
plt.title('Comparison at point (' + str(point_x) + ',' + str(point_y) + ')')
plt.grid()
plt.savefig(inputdir + save_name[index] + '_comparison_for2_adh_sswm_' + str(point_x) + '_' + str(point_y) + '.pdf', bbox_inches = 'tight')
plt.close()


sswm, adh = [], []

with open(inputdir + file[0], 'r') as f:
    for j, line in enumerate(f):
        if j != 0:
            line = line.split(',')
            sswm.append(float(line[0]))

with open(inputdir + file[1], 'r') as f:
    for j, line in enumerate(f):
        if j != 0:
            line = line.split(',')
            adh.append(float(line[0]))

for i in [sswm[j]-adh[j] for j in range(len(sswm))]:
    print(i)

print(max([sswm[j]-adh[j] for j in range(len(sswm))]))
