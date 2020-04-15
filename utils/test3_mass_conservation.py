import numpy as np
import matplotlib
matplotlib.use('Agg')

font = {'family': 'normal', 'size': 22}
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

out_put_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/test3_res/SSWM_results/"
mass_file = "total_mass_at_every_time_step.npy"
time_file = "time_stamp_at_every_time_step.npy"
scale_max = 1.05
scale_min = 0.95
scale_in_max = 1.0005
scale_in_min = 0.9995
num = 6
marker = 2
steptime = 1

mass = np.load(out_put_dir + mass_file)
time = np.load(out_put_dir + time_file)

mass = mass[::steptime]
time = time[::steptime]

print(mass)
print(time)

#pdb.set_trace()
plt.figure(figsize=[17, 7])
plt.plot(time, mass, '-o', markersize=marker)
plt.xlim([min(time), max(time)])
plt.xticks(np.linspace(min(time), max(time), num))
plt.ylim([min(mass)*scale_in_min, max(mass)*scale_in_max])
plt.yticks(np.linspace(min(mass)*scale_in_min, max(mass)*scale_in_max, num))
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
plt.grid()
plt.xlabel("time: sec")
plt.ylabel("volume of mass: " + r'$m^3$')
plt.title("Conservation of mass")
plt.legend(["zoom in"])
# plt.show()
plt.savefig(out_put_dir + "Conservation_of_mass_zoom_in_6_hr.pdf", bbox_inches='tight')

plt.figure(figsize=[17, 7])
plt.plot(time, mass, '-o', markersize=marker)
plt.ylim([min(mass)*scale_min, max(mass)*scale_max])
plt.yticks(np.linspace(min(mass)*scale_min, max(mass)*scale_max, num))
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
plt.grid()
plt.xlim([min(time), max(time)])
plt.xticks(np.linspace(min(time), max(time), num))
plt.xlabel("time: sec")
plt.ylabel("volume of mass: " + r'$m^3$')
plt.title("Conservation of mass")
plt.legend(["zoom out"])
# plt.show()
plt.savefig(out_put_dir + "Conservation_of_mass_zoom_out_6_hr.pdf", bbox_inches='tight')
