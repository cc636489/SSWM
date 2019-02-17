import numpy as np
import sys
sys.path.append('/workspace/Documentation/Research_Doc/SFEM_Doc/7-NS-bitbucket/src/')
from makestobasis import *
import chaospy as cp
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mtick


loc = '/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/'
dir = 'test2_res/'
name = 'u'
file = loc + dir + name + '_modes_25_50_ts.csv'
scale = 1.05
ns = 1000


distname = "uniform"
sto_poly_deg = 4
sto_poly_dim = 2
coefficient = [-1.2, 1.2, -2, 2]

basis_struct = MakeStoBasis(distname, sto_poly_deg, sto_poly_dim, coefficient)
orthpol = basis_struct.get("basis")
JointCDF = basis_struct.get("jointcdf")
nmodes = basis_struct.get('nmodes')

# check orthogonality
for i in range(len(orthpol)):
    print cp.E(orthpol[i]*orthpol[i], JointCDF)

# read in row, col
with open(file, 'r') as f:
    for i, line in enumerate(f):
        line = line.split(',')
        col = len(line) - 1
    row = i

# initialization
time = np.arange(0, row, 1)
index = np.zeros(col, dtype=int)
field = np.zeros([row, col])
field_std = np.zeros(row)

# file read
with open(file, 'r') as f:
    for i, line in enumerate(f):
        line = line.split(',')
        if i == 0:
            for j in range(col):
                index[j] = int(line[j+1])
        if i != 0:
            # time[i-1] = float(line[0])
            for j in range(col):
                field[i-1, j] = float(line[j+1])

# model formulate
sample = JointCDF.sample(ns)

for step in range(1, row):
    random_res = []
    for j in range(ns):
        sum = 0
        for l, i in enumerate(index):
            sum += field[step, l] * orthpol[i](sample[0, j], sample[1, j])
        random_res.append(sum)
    # plot pdf
    plt.figure(figsize=[7, 6])
    sns.distplot(random_res, bins=30)
    plt.xlabel(name)
    plt.ylabel('pdf')
    plt.xticks(np.linspace(min(random_res), max(random_res), 8))
    plt.gca().xaxis.set_major_formatter(mtick.FormatStrFormatter('%.3f'))
    plt.title('time = ' + str(step) + ' sec')
    # plt.show()
    plt.savefig(file[:-4]+'_'+str(step)+'.png')
