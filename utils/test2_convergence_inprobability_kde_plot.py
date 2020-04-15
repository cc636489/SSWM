
import numpy as np
import matplotlib
matplotlib.use('Agg')

font = {'family' : 'normal',
              'size'   : 22}

matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mtick


input_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/test2_res/convergence_inprobability_test_result_dimension_2/"
name = "u"
point_x = 25.0
point_y = 25.0
color = ['mediumseagreen', 'darkorange', 'deepskyblue', 'darkred',  'orchid', 'slategray', 'black']
order_list = [0,1,2,3]

plt.subplots(figsize=(11.5, 7))

for i,order in enumerate(order_list):
    file_name = "binrandom_surrogate_"+name.lower()+"_samples_order_"+str(point_x)+"_"+str(point_y)+"_"+str(order)+".npy"
    ensembles = np.load(input_dir+file_name)
    sns.kdeplot(ensembles,label="order="+str(order), color=color[i],linestyle='-')
    plt.axvline(x=np.mean(ensembles), color=color[i], linewidth = 0.5, linestyle='dashed')
    plt.legend(loc="upper right")

plt.savefig(input_dir+name+"_convergence_inprobability_plot_"+str(point_x)+"_"+str(point_y)+".pdf")
