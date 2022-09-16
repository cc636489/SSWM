
import numpy as np
import matplotlib
matplotlib.use('Agg')

font = {'family' : 'normal',
              'size'   : 22}

matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mtick


input_dir = "/Users/chenchen/gloria/test_fenics/SupportingRun/test2_convergence_at_space_25.0_25.0_time_1s_dimension_1/"
name = "eta"
point_x = 25.0
point_y = 25.0
color = ['mediumseagreen', 'darkorange', 'deepskyblue', 'darkred',  'orchid', 'slategray', 'black']
order_list = [1,2,3]

plt.subplots(figsize=(11.5, 7))

for i,order in enumerate(order_list):
    file_name = "binrandom_surrogate_"+name.lower()+"_samples_order_"+str(point_x)+"_"+str(point_y)+"_"+str(order)+".npy"
    ensembles = np.load(input_dir+file_name)
#    sns.kdeplot(ensembles,label="order="+str(order), color=color[i],linestyle='-')
    sns.distplot(ensembles,label="order="+str(order), color=color[i])
    plt.axvline(x=np.mean(ensembles), color=color[i], linewidth = 0.5, linestyle='dashed')
    plt.legend(loc="upper right")

plt.savefig(input_dir+name+"_convergence_inprobability_plot_"+str(point_x)+"_"+str(point_y)+".pdf")
