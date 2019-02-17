import numpy as np
import matplotlib.pyplot as plt

loc = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/'
dir = 'test2_res/'
name = 'u'
file = loc + dir + name + '_modes_25_50_ts.csv'
scale = 1.05

with open(file, 'r') as f:
    for i, line in enumerate(f):
        line = line.split(',')
        col = len(line) - 2
    row = i

time = np.arange(0, row, 1)
mean = np.zeros(row)
field = np.zeros([row, col])
field_std = np.zeros(row)

with open(file, 'r') as f:
    for i, line in enumerate(f):
        line = line.split(',')
        if i != 0:
            # time[i-1] = float(line[0])
            mean[i-1] = float(line[1])
            for j in range(col):
                field[i-1, j] = float(line[j+2])

for i in range(row):
    sum = 0
    for j in range(col):
        sum += field[i, j] ** 2
    field_std[i] = np.sqrt(sum)

plt.figure(figsize=[7, 6])
plt.plot(time, mean, '-*', alpha=0.2, color='#1B2ACC')
plt.fill_between(time, mean-field_std, mean+field_std, alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF',
    linewidth=0.2, linestyle='dashdot', antialiased=True)
# plt.errorbar(time, mean, yerr=field_std, fmt='-^', color='#FF9848')
plt.xlim([min(time), max(time)*scale])
plt.ylim([min(mean-field_std), max(mean+field_std)*scale])
plt.xticks(np.linspace(min(time), max(time)*scale, 11))
plt.yticks(np.linspace(min(mean-field_std), max(mean+field_std)*scale, 11))
plt.xlabel('time(sec)')
plt.ylabel(name)
plt.legend(['mean', r'mean$\pm\sigma$'])
# plt.show()
plt.savefig(file[:-4]+'.png')
