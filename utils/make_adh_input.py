import numpy as np

A = 0.75
period = 12.41664*60*60
omega = 2*np.pi/period
time = np.linspace(0, 223499.52000, 5000)
ele = np.zeros(len(time))
for i,t in enumerate(time):
    ele[i] = A*np.sin(omega*t)
    print time[i],"  ", ele[i]

