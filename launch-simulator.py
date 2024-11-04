import numpy as np
import matplotlib.pyplot as plt
import atmospheric_models
import rocket_models
import general_functions

# ---------------------ASSUMPTIONS-------------------- #
# 1. Earth is a perfect sphere with uniform density (point-mass gravity model)
# 2. A day is exactly 24 hours
# 3. The Rocket is radially symmetrical (all force locations on centerline)
# 4. Negligible wind, all particles of air move in a perfect circle about the earth's axis of rotation each day
# 5. Center of Pressure is static during each stage
# 6. Liquid Fuel has no slosh, and drains in the form of a perfect cylinder of decreasing height
# 7. Ground effect is not considered
# 8. Lift force is negligible (C_l = 0.0)
# 9. Mass, dimension properties and MOIs are simplified as defined in rocket_models.py


# Environment Parameters
dt = 1.0                            # seconds

# Earth Parameters
G = 6.67430*10**-11                 # Earth Parameters - Point Mass Model
r_earth = 6378137                   # meters
m_earth = 5.9722e24                 # kg
thetaDot_earth = 2*np.pi/60/60/24   # radians/sec earth rotation

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

#Establish Models
atmospheric_model = atmospheric_models.US_Standard_Atmosphere()
vehicle_model = rocket_models.FalconIX()

# Vehicle Parameters
velocity = [8000.0, 0, 0.0]
coords = [0.0, r_earth+250000.0, 0.0]
twist = [0.0, 0.0, 0.0]
# vehicle_quaternion = [0.5, 0.5, 0.5, 0.5] # directly in +Y direction TODO: Convert all rotations to quaternion
body_twist = [0.0, 0.0, 0.0]
drag_cross_section = vehicle_model.cross_section_area

# Instantiate Debugging Data Arrays
t = []

Ftx_log = []
Fty_log = []
Ftz_log = []

Fdx_log = []
Fdy_log = []
Fdz_log = []

Fgx_log = []
Fgy_log = []
Fgz_log = []

x_log = []
y_log = []
z_log = []

vx_log = []
vy_log = []
vz_log = []

ax_log = []
ay_log = []
az_log = []

mass_log = []


# Determine Orbital Period
semimajor_axis = (2/np.linalg.norm(coords)-np.linalg.norm(velocity)**2/(G*m_earth))**-1
orbital_period = 2*np.pi*np.sqrt(semimajor_axis**3/(G*m_earth))

# for i in np.linspace(0,orbital_period,int(orbital_period/dt)):
for i in np.linspace(0,2000,int(2000/dt)):
    dist = np.linalg.norm(coords)
    if dist < r_earth:
        print('Rocket has landed rapidly (crashed)')
        break

    pressure,temperature,density = atmospheric_model.get_params(dist)


    F_thrust, body_thrust_vector, engine_distance = vehicle_model.getThrustVector(twist, pressure)
    thrust_vector = general_functions.inverseCoordinateTransform(body_thrust_vector, twist) # TODO: Verify This
    print(thrust_vector,'.....................',body_thrust_vector)
    # Drag Calculations
    air_ground_velocity = np.cross([0,0,thetaDot_earth],[coords[0],coords[1],0])
    air_speed_vector = velocity + air_ground_velocity
    air_speed = np.linalg.norm(air_speed_vector)
    Cd = vehicle_model.getDragCoefficient(air_speed,density,pressure)
    F_drag = 0.5*density*drag_cross_section*Cd*air_speed**2
    
    # Gravity Calculations
    vehicle_mass = vehicle_model.getMass()
    F_gravity = G*m_earth*vehicle_mass/(dist**2)
    gravity_vector = np.negative(coords)/dist

    # Update Step
    Fg_x = gravity_vector[0]*F_gravity
    Fg_y = gravity_vector[1]*F_gravity
    Fg_z = gravity_vector[2]*F_gravity

    Fd_x = -F_drag*air_speed_vector[0]
    Fd_y = -F_drag*air_speed_vector[1]
    Fd_z = -F_drag*air_speed_vector[2]

    Ft_x = F_thrust*thrust_vector[0]
    Ft_y = F_thrust*thrust_vector[1]
    Ft_z = F_thrust*thrust_vector[2]

    forces_x = Fg_x + Fd_x + Ft_x
    forces_y = Fg_y + Fd_y + Ft_y
    forces_z = Fg_z + Fd_z + Ft_z
    
    accel_x = forces_x/vehicle_mass*dt
    accel_y = forces_y/vehicle_mass*dt
    accel_z = forces_z/vehicle_mass*dt

    velocity[0] += accel_x
    velocity[1] += accel_y
    velocity[2] += accel_z

    # Transform Forces to Body Frame
    Fg_rotated = general_functions.coordinateTransform([Fg_x,Fg_y,Fg_z],twist)
    Fd_rotated = general_functions.coordinateTransform([Fd_x,Fd_y,Fd_z],twist)
    Ft_rotated = general_functions.coordinateTransform([Ft_x,Ft_y,Ft_z],twist)

    # Get Distances of Forces with respect to payload
    center_of_gravity = vehicle_model.getGravityCenter()
    center_of_pressure = vehicle_model.getCenterofPressure()

    Ixx, Iyy, Izz = vehicle_model.getMOIs()

    moments_xx = Fg_rotated[0]*center_of_gravity + Fd_rotated[0]*center_of_pressure + Ft_rotated[0]*engine_distance
    moments_yy = Fg_rotated[1]*center_of_gravity + Fd_rotated[1]*center_of_pressure + Ft_rotated[1]*engine_distance
    moments_zz = Fg_rotated[2]*center_of_gravity + Fd_rotated[2]*center_of_pressure + Ft_rotated[2]*engine_distance

    # Inertial Frame Rotations and Coordinates
    body_twist[0] = moments_xx/Ixx*dt  
    body_twist[1] = moments_yy/Iyy*dt
    body_twist[2] = moments_zz/Izz*dt

    # twist += general_functions.angleRates_BodytoEuler(body_twist,twist)  # TODO: Verify this, highly suspect

    coords[0] += velocity[0]*dt 
    coords[1] += velocity[1]*dt 
    coords[2] += velocity[2]*dt

    # Data Logging
    t.append(i)
    x_log.append(coords[0]) 
    y_log.append(coords[1]) 
    z_log.append(coords[2])

    Ftx_log.append(Ft_x)
    Fty_log.append(Ft_y)
    Ftz_log.append(Ft_z)

    Fdx_log.append(Fd_x)
    Fdy_log.append(Fd_y)
    Fdz_log.append(Fd_z)

    Fgx_log.append(Fg_x)
    Fgy_log.append(Fg_y)
    Fgz_log.append(Fg_z)

    vx_log.append(velocity[0])
    vy_log.append(velocity[1])
    vz_log.append(velocity[2])

    ax_log.append(accel_x)
    ay_log.append(accel_y)
    az_log.append(accel_z)

    mass_log.append(vehicle_mass)


# Configure 3-D Plot
frame_multiplier = 1.5
ax.set_xlim([-frame_multiplier*r_earth,frame_multiplier*r_earth])
ax.set_ylim([-frame_multiplier*r_earth,frame_multiplier*r_earth])
ax.set_zlim([-frame_multiplier*r_earth,frame_multiplier*r_earth])
ax.set_box_aspect([1,1,1])
ax.plot(x_log,y_log,z_log)
ax.set_xlabel('X [m]')
ax.set_ylabel('Y [m]')
ax.set_zlabel('Z [m]')

# Generate Translucent Earth
u, v = np.mgrid[0:2*np.pi:64j, 0:np.pi:64j]
x_earth = r_earth*np.cos(u)*np.sin(v)
y_earth = r_earth*np.sin(u)*np.sin(v)
z_earth = r_earth*np.cos(v)
ax.plot_surface(x_earth, y_earth, z_earth, color="b",alpha=0.25)
# plt.show()

# Generate Debugging Plots
fig2, ax2 = plt.subplots(nrows = 3, ncols = 3)
ax2[0,0].plot(t,Ftx_log)
ax2[0,0].set_title('Thrust Forces')
ax2[0,0].set_xlabel('time [s]')
ax2[0,0].set_ylabel('Thrust Force X [N]')
ax2[1,0].plot(t,Fty_log)
ax2[1,0].set_xlabel('time [s]')
ax2[1,0].set_ylabel('Thrust Force Y [N]')
ax2[2,0].plot(t,Ftz_log)
ax2[2,0].set_xlabel('time [s]')
ax2[2,0].set_ylabel('Thrust Force Z [N]')
ax2[0,1].plot(t,Fdx_log)
ax2[0,1].set_title('Drag Forces Forces')
ax2[0,1].set_xlabel('time [s]')
ax2[0,1].set_ylabel('Drag Force X [N]')
ax2[1,1].plot(t,Fdy_log)
ax2[1,1].set_xlabel('time [s]')
ax2[1,1].set_ylabel('Drag Force Y [N]')
ax2[2,1].plot(t,Fdz_log)
ax2[2,1].set_xlabel('time [s]')
ax2[2,1].set_ylabel('Drag Force Z [N]')
ax2[0,2].plot(t,Fgx_log)
ax2[0,2].set_title('Gravitational Forces')
ax2[0,2].set_xlabel('time [s]')
ax2[0,2].set_ylabel('Gravitational Force X [N]')
ax2[1,2].plot(t,Fgy_log)
ax2[1,2].set_xlabel('time [s]')
ax2[1,2].set_ylabel('Gravitational Force Y [N]')
ax2[2,2].plot(t,Fgz_log)
ax2[2,2].set_xlabel('time [s]')
ax2[2,2].set_ylabel('Gravitational Force Z [N]')
plt.tight_layout()

# Trajectory Debugging Plots
fig3, ax3 = plt.subplots(nrows = 3, ncols = 3)
ax3[0,0].plot(t,x_log)
ax3[0,0].set_title('Positions')
ax3[0,0].set_xlabel('time [s]')
ax3[0,0].set_ylabel('X [m]')
ax3[1,0].plot(t,y_log)
ax3[1,0].set_xlabel('time [s]')
ax3[1,0].set_ylabel('Y [m]')
ax3[2,0].plot(t,z_log)
ax3[2,0].set_xlabel('time [s]')
ax3[2,0].set_ylabel('Z [m]')
ax3[0,1].plot(t,vx_log)
ax3[0,1].set_title('Velocities')
ax3[0,1].set_xlabel('time [s]')
ax3[0,1].set_ylabel('Velocity X [m/s]')
ax3[1,1].plot(t,vy_log)
ax3[1,1].set_xlabel('time [s]')
ax3[1,1].set_ylabel('Velocity Y [m/s]')
ax3[2,1].plot(t,vz_log)
ax3[2,1].set_xlabel('time [s]')
ax3[2,1].set_ylabel('Velocity Z [m]')
ax3[0,2].plot(t,ax_log)
ax3[0,2].set_title('Accelerations')
ax3[0,2].set_xlabel('time [s]')
ax3[0,2].set_ylabel('Acceleration X [m/s]')
ax3[1,2].plot(t,ay_log)
ax3[1,2].set_xlabel('time [s]')
ax3[1,2].set_ylabel('Acceleration Y [m/s]')
ax3[2,2].plot(t,az_log)
ax3[2,2].set_xlabel('time [s]')
ax3[2,2].set_ylabel('Acceleration Z [m]')
plt.tight_layout()

fig4, ax4 = plt.subplots(nrows = 2, ncols = 1)
ax4[0].plot(t,mass_log)
plt.show()