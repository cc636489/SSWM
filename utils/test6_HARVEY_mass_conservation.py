import numpy as np
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.patches import ConnectionPatch
matplotlib.use('Agg')

font = {'family': 'normal', 'size': 22}
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick

out_put_dir = "/Users/chenchen/gloria/4-NS-results-and-tests/regression_test/SSWM_comparison_Gulf_wind_IKE/IKE_correct_bug_with_atmos_press/"
# out_put_dir = "/workspace/Documentation/Research_Doc/SFEM_Doc/4-NS-results-and-tests/regression_test/SSWM_comparison_Gulf_wind_HARVEY/HARVEY_correct_bug_with_atmos_press/"
mass_file = "total_mass_at_every_time_step.npy"
time_file = "time_stamp_at_every_time_step.npy"
scale_max = 1.01
scale_min = 0.95
num = 6
marker = 8
steptime = 10 

mass = np.load(out_put_dir + mass_file)
time = np.load(out_put_dir + time_file)

mass = mass[::steptime]
time = time[::steptime]/24.0/3600.0

# mass = mass[30:45:1]
# time = time[30:45:1]

print(mass)
print(time)

#pdb.set_trace()
# plt.figure(figsize=[10, 7])
# plt.plot(time, mass, '-*', markersize=marker, color='mediumseagreen')
# plt.xlim([min(time), max(time)])
# plt.xticks(np.linspace(min(time), max(time), num))
# plt.ylim([min(mass)-10, max(mass)+10])
# # plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
# plt.xlabel("time: days")
# plt.ylabel("volume of mass: " + r'$m^3$')
# # plt.title("Conservation of mass")
# plt.legend(["zoom in"])
# # plt.show()
# plt.savefig(out_put_dir + "Conservation_of_mass_zoom_in_6_hr_latest.pdf")

# plt.figure(figsize=[17, 7])
# plt.plot(time, mass, '-*', markersize=marker, color='mediumseagreen')
# plt.ylim([min(mass)*scale_min, max(mass)*scale_max])
# plt.yticks(np.linspace(min(mass)*scale_min, max(mass)*scale_max, num))
# plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
# plt.xlim([min(time), max(time)])
# plt.xticks(np.linspace(min(time), max(time), num))
# plt.xlabel("time: days")
# plt.ylabel("volume of mass: " + r'$m^3$')
# plt.title("Conservation of mass")
# plt.legend(["zoom out"])
# # plt.show()
# plt.savefig(out_put_dir + "Conservation_of_mass_zoom_out_6_hr_include_small_pic.pdf")


# 画局部放大图

fig, ax = plt.subplots(1, 1, figsize=[17, 7])
ax.plot(time, mass, '-*', markersize=marker, color='mediumseagreen')
ax.set_ylim([min(mass)*scale_min, max(mass)*scale_max])
ax.set_yticks(np.linspace(min(mass)*scale_min, max(mass)*scale_max, num))
plt.gca().yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
ax.set_xlim([min(time), max(time)])
ax.set_xticks(np.linspace(min(time), max(time), num))
ax.set_xlabel("time: days")
ax.set_ylabel("volume of mass: " + r'$m^3$')
ax.set_title("Conservation of mass")
ax.legend(["zoom out"], loc="lower right")
# plt.show()
fig.savefig(out_put_dir + "Conservation_of_mass_zoom_out_6_hr_latest.pdf")

axins = ax.inset_axes((0.2, 0.2, 0.4, 0.3))

