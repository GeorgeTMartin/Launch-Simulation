import time, os
import numpy as np
from MPC_trajectory_optimization_code_generator import *
import matplotlib.pyplot as plt
'''
STATES: x_dot, y_dot, x_d_dot, y_d_dot, theta, t_dot
INPUTS: theta
'''
G = 6.67430*10**-11                                             # Earth Parameters - Point Mass Model
earth_radius = 6378137                                          # meters
earth_mass = 5.9722e24                                          # kg
orbit_altitude = 0.5e9                                          # meters
orbit_speed = np.sqrt(G*earth_mass/orbit_altitude)              # m/s 
Tf = 600  # prediction horizon
N = 600  # number of discretization steps

# Load Model
constraint, model, acados_solver = acados_generator(rebuild_flag=True)

# Dimensions
nx = 6 # model.x.size()[0]
nu = 1 # model.u.size()[0]
# print("size nx, nu: {} | {}".format(nx,nu))
ny = nx + nu

Nsim = N

# initialize data structs
simX = np.ndarray((Nsim, nx))
simU = np.ndarray((Nsim, nu))
tcomp_sum = 0
tcomp_max = 0
x0 = model.x0
acados_solver.set(0, "lbx", x0)
acados_solver.set(0, "ubx", x0)

# simulate
# update reference
for j in range(N):
    try:
        yref = np.array([0, orbit_altitude, orbit_speed, 0, 0, 0, 0])
        acados_solver.set(j, "yref", yref)
    except:
        yref = np.array([0, orbit_altitude, orbit_speed, 0, 0, 0])
        acados_solver.set(j, "yref", yref)
    

        
yref_N = yref[:-1]
acados_solver.set(N, "yref", yref_N)

# solve ocp
acados_solver.solve()

# get solution
x = []
y = []
x_dot = []
y_dot = []
theta = []
t = []
u = []

for i in range(N):
    temp_state = acados_solver.get(i, "x")
    x.append(temp_state[0])
    y.append(temp_state[1])
    x_dot.append(temp_state[2])
    y_dot.append(temp_state[3])
    theta.append(temp_state[4])
    u.append(temp_state[5])

t=np.linspace(0,Tf,N)
fig, ax = plt.subplots(6,1)
ax[0].plot(t,x)
ax[1].plot(t,y)
ax[2].plot(t,x_dot)
ax[3].plot(t,y_dot)
ax[4].plot(t,theta)
ax[5].plot(t,u)
plt.show()
