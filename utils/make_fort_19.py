import numpy as np

# Inlet
# tidal_amplitude = 0.75  # galveston bay area normal tides.
# tidal_period = 12.41664*60*60
# start_time = 0.0
# end_time = tidal_period * 5
# time_step = tidal_period/100
# omega = 2*np.pi/tidal_period
# time = np.linspace(0, 223500, 501)


# hump
NETA = 21
Time_Increment = 1
tidal_amplitude = 0.1
omega = np.pi/20.
time = np.linspace(0, 401, 402)
ele = [0]*len(time)

f2 = open('/workspace/RMS_runs/Hump test case/fort.19_sin_func', "w+")
print Time_Increment
f2.write(str(Time_Increment))
f2.write('\n')
for i, t in enumerate(time):

    ele[i] = tidal_amplitude*np.sin(omega*t)
    for _ in range(NETA):
        print ele[i]
        f2.write(str(ele[i]))
        f2.write('\n')

f2.close()

print time