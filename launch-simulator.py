import numpy as np
import matplotlib.pyplot as plt
import atmospheric_models
import rocket_models

G = 6.67430*10**-11# Earth Parameters - Point Mass Model
r_earth = 6378137 #meters
m_earth = 5.9722e24 #kg
thetaDot_earth = 2*np.pi/60/24 #radians/sec earth rotation
# Satellite Parameters
velocity = [7500/np.sqrt(2), 0, 7500/np.sqrt(2)]
coords = [0,r_earth+1000000,0] # Initial State
m_satellite = 100 #kg
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
x0_plot = []
y0_plot = []
z0_plot = []

x_plot = []
y_plot = []
z_plot = []
dt = 0.1 # seconds

#Determine Orbital Period
mu = np.linalg.norm(velocity)**2*np.linalg.norm(coords)
semimajor_axis =(2/np.linalg.norm(coords)-np.linalg.norm(velocity)**2/(G*m_earth))**-1
orbital_period = 2*np.pi*np.sqrt(semimajor_axis**3/(G*m_earth))

#Establish Models
atmospheric_model = atmospheric_models.US_Standard_Atmosphere()
rocket_model = rocket_models.LiquidRocket()

for i in np.linspace(0,orbital_period,int(orbital_period/dt)):
    dist = np.linalg.norm(coords)
    if dist < r_earth:
        print('Rocket has landed rapidly (crashed)')
        break

    [pressure,temperature, density] = atmospheric_model.get_params(dist)

    # Determine Drag Coefficient
   
    air_relative_speed = 0.0
    air_speed = 0.0
    Cd = rocket_model.getDragCoefficient(air_speed,density,pressure)

    F_gravity = G*m_earth*m_satellite/(dist**2)
    gravity_vector = np.negative(coords)/dist
    velocity[0] += gravity_vector[0]*F_gravity/m_satellite*dt
    velocity[1] += gravity_vector[1]*F_gravity/m_satellite*dt
    velocity[2] += gravity_vector[2]*F_gravity/m_satellite*dt

    coords[0] += velocity[0]*dt 
    coords[1] += velocity[1]*dt 
    coords[2] += velocity[2]*dt

    x_plot.append(coords[0]) 
    y_plot.append(coords[1]) 
    z_plot.append(coords[2])



frame_multiplier = 1.5
ax.set_xlim([-frame_multiplier*r_earth,frame_multiplier*r_earth])
ax.set_ylim([-frame_multiplier*r_earth,frame_multiplier*r_earth])
ax.set_zlim([-frame_multiplier*r_earth,frame_multiplier*r_earth])
ax.set_box_aspect([1,1,1])
ax.plot(x_plot,y_plot,z_plot)
ax.plot(x0_plot,y0_plot,z0_plot)
u, v = np.mgrid[0:2*np.pi:64j, 0:np.pi:64j]
x_earth = r_earth*np.cos(u)*np.sin(v)
y_earth = r_earth*np.sin(u)*np.sin(v)
z_earth = r_earth*np.cos(v)
ax.plot_surface(x_earth, y_earth, z_earth, color="b",alpha=0.25)
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_zlabel('Z [m]')
plt.show()
