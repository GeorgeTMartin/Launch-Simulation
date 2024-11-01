import numpy as np
import matplotlib.pyplot as plt
import atmospheric_models
import rocket_models

# ---------------------ASSUMPTIONS-------------------- #
# 1. Earth is a perfect sphere with uniform density (point-mass gravity model)
# 2. A day is exactly 24 hours
# 3. The Rocket is radially symmetrical (all force locations on centerline)
# 4. Negligible wind, all particles of air move in a perfect circle about the earth's axis of rotation each day
# 5. Center of Pressure is static during each stage
# 6. Liquid Fuel has no slosh, and drains in the form of a perfect cylinder of decreasing height
# 7. Ground effect is not considered
# 8. Lift force is negligible


# Environment Parameters
dt = 1.0                            # seconds

# Earth Parameters
G = 6.67430*10**-11                 # Earth Parameters - Point Mass Model
r_earth = 6378137                   # meters
m_earth = 5.9722e24                 # kg
thetaDot_earth = 2*np.pi/60/60/24   # radians/sec earth rotation

fig = plt.figure()
ax = fig.add_subplot(projection='3d')
x0_plot = []
y0_plot = []
z0_plot = []

x_plot = []
y_plot = []
z_plot = []

#Establish Models
atmospheric_model = atmospheric_models.US_Standard_Atmosphere()
vehicle_model = rocket_models.OrbitingSatellite()

# Vehicle Parameters
velocity = [7500.0, 0.0, 0.0]
coords = [0.0, r_earth+1000000.0, 0.0]
twist = [0.0, 0.0, 0.0]

# Instantiate Debugging Data Arrays
t = []
testpointx = []
testpointy = []

# Determine Orbital Period
semimajor_axis = (2/np.linalg.norm(coords)-np.linalg.norm(velocity)**2/(G*m_earth))**-1
orbital_period = 2*np.pi*np.sqrt(semimajor_axis**3/(G*m_earth))

for i in np.linspace(0,orbital_period,int(orbital_period/dt)):
    dist = np.linalg.norm(coords)
    if dist < r_earth:
        print('Rocket has landed rapidly (crashed)')
        break

    pressure,temperature,density = atmospheric_model.get_params(dist)


    F_thrust, thrust_vector, engine_distance = vehicle_model.getThrustVector(twist)
    
    # Drag Calculations
    air_ground_velocity = np.cross([0,0,thetaDot_earth],[coords[0],coords[1],0])
    air_speed_vector = velocity + air_ground_velocity
    air_speed = np.linalg.norm(air_speed_vector)
    Cd = vehicle_model.getDragCoefficient(air_speed,density,pressure)
    F_drag = 0.5*density*Cd*air_speed**2
    
    # Gravity Calculations
    vehicle_mass = vehicle_model.getMass()
    F_gravity = G*m_earth*vehicle_mass/(dist**2)
    gravity_vector = np.negative(coords)/dist

    # Update Step
    Fg_x = gravity_vector[0]*F_gravity
    Fg_y = gravity_vector[1]*F_gravity
    Fg_z = gravity_vector[2]*F_gravity

    Fd_x = -F_drag*velocity[0]
    Fd_y = -F_drag*velocity[1]
    Fd_z = -F_drag*velocity[2]

    Ft_x = F_thrust*thrust_vector[0]
    Ft_y = F_thrust*thrust_vector[1]
    Ft_z = F_thrust*thrust_vector[2]

    forces_x = Fg_x + Fd_x + Ft_x
    forces_y = Fg_y + Fd_y + Ft_y
    forces_z = Fg_z + Fd_z + Ft_z
    
    velocity[0] += forces_x/vehicle_mass*dt
    velocity[1] += forces_y/vehicle_mass*dt
    velocity[2] += forces_z/vehicle_mass*dt

    # Moment Calculations TODO: Add Moment Calculations NOTE: Torque = radial vector cross force vector
    center_of_gravity = vehicle_model.getGravityCenter()

    moments_xx = 0.0
    moments_yy = 0.0
    moments_zz = 0.0

    Ixx, Iyy, Izz = vehicle_model.getMOIs()

    rotation_x_rate = moments_zz/Izz*dt
    rotation_y_rate = moments_zz/Iyy*dt
    rotation_z_rate = moments_zz/Izz*dt

    # Inertial Frame Rotations and Coordinates
    twist[0] = 0.0
    twist[1] = 0.0
    twist[2] = 0.0

    coords[0] += velocity[0]*dt 
    coords[1] += velocity[1]*dt 
    coords[2] += velocity[2]*dt

    # Data Logging
    t.append(i)
    x_plot.append(coords[0]) 
    y_plot.append(coords[1]) 
    z_plot.append(coords[2])

    testpointx.append(Cd)
    testpointy.append(coords[0])



# Configure 3-D Plot
frame_multiplier = 1.5
ax.set_xlim([-frame_multiplier*r_earth,frame_multiplier*r_earth])
ax.set_ylim([-frame_multiplier*r_earth,frame_multiplier*r_earth])
ax.set_zlim([-frame_multiplier*r_earth,frame_multiplier*r_earth])
ax.set_box_aspect([1,1,1])
ax.plot(x_plot,y_plot,z_plot)
ax.plot(x0_plot,y0_plot,z0_plot)
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_zlabel('Z [m]')

# Generate Translucent Earth
u, v = np.mgrid[0:2*np.pi:64j, 0:np.pi:64j]
x_earth = r_earth*np.cos(u)*np.sin(v)
y_earth = r_earth*np.sin(u)*np.sin(v)
z_earth = r_earth*np.cos(v)
ax.plot_surface(x_earth, y_earth, z_earth, color="b",alpha=0.25)
plt.show()

# Generate Debugging Plots
fig2 = plt.figure()
ax2 = fig2.add_subplot()
ax3 = ax2.twinx()
ax2.plot(t,testpointx,'#FFAE42')
ax2.set_xlabel('Time [s]')
ax3.plot(t,testpointy,'#FF5349')
ax2.set_ylabel('Test Point X',color='#FFAE42')
ax3.set_ylabel('Test Point Y',color='#FF5349')
plt.show()