
import numpy as np
import matplotlib
matplotlib.use('Agg')

font = {'family' : 'normal',
              'size'   : 22}

matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.ticker as mtick


input_dir = "/Users/chenchen/gloria/test_fenics/stochastic/test4_stochastic_inlet/"
truth_dir = "/Users/chenchen/gloria/4-NS-results-and-tests/regression_test_stochastic/test4_bins/"
name = "eta"
xlabel = "Surface Elevation(m)"
#name = "u"
#xlabel = "x-direction velocity(m/s)"
point_x = 0.0
point_y = 0.0
timestep = 405
color = ['mediumseagreen', 'darkorange', 'deepskyblue', 'darkred',  'orchid', 'slategray', 'black']
order_list = [3]

plt.subplots(figsize=(11.5, 7))

for i,order in enumerate(order_list):
    file_name = "binrandom_surrogate_"+name.lower()+"_samples_order_"+str(point_x)+"_"+str(point_y)+"_"+str(timestep)+"_NSample_5000.npy"
    truth_name = "test4_bin_random_"+name+"1_true_all_points_order_1.npy"
    ensembles = np.load(input_dir+file_name)
    truth = np.load(truth_dir+truth_name)

    #sns.kdeplot(ensembles,label="order="+str(order), color=color[i],linestyle='-')
    sns.distplot(ensembles,label="SSWM", color=color[i])
    sns.distplot(truth[:, timestep, 1],label="truth", color=color[i+1]) # timestep=40, testnode number = 0.
    ensemble_mean = np.mean(ensembles)
    ensemble_var = np.sqrt(np.var(ensembles))
    truth_mean = np.mean(truth[:, timestep, 1])
    truth_var = np.sqrt(np.var(truth[:, timestep, 1]))
    print("test:============= "+name+" ========"+str(timestep)+"==========")
    print("ensemble mean", ensemble_mean)
    print("truth mean", truth_mean)
    print("ensemble devia", ensemble_var)
    print("truth devia", truth_var)
#    plt.axvline(x=ensemble_mean, label="SSWM mean", color=color[i], linewidth = 0.8, linestyle='dashed')
#    plt.axvline(x=truth_mean, label="truth mean", color=color[i+1], linewidth = 0.8, linestyle='dashed')
    plt.plot([ensemble_mean - ensemble_var, ensemble_mean, ensemble_mean + ensemble_var],[0,0,0], color=color[i], marker="o",markersize=22, label = "SSWM:"+r'$\mu\pm\sigma$')
    plt.plot([truth_mean - truth_var, truth_mean, truth_mean + truth_var],[0,0,0], marker="*", color=color[i+1], markersize=22, label = "truth:"+r'$\mu\pm\sigma$')
    plt.xlabel(xlabel)
    #plt.title("time:1.061 days")
    plt.title("time:2.095 days")
    plt.legend(loc="center")

plt.savefig(input_dir+name+"_pdf_comparison_plot_"+str(point_x)+"_"+str(point_y)+"_timeStep_"+str(timestep)+".png")
