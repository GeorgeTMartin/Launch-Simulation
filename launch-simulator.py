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
# 10. Fairing inertial properties not considered

# Simulation Environment Variables
dt = 1.0 # Step Time [s]
simulation_time = 10000 # [s]
atmospheric_model = atmospheric_models.US_Standard_Atmosphere()
vehicle_model = rocket_models.FalconIX()

# Earth Parameters
G = 6.67430*10**-11                                             # Earth Parameters - Point Mass Model
earth_radius = 6378137                                          # meters
earth_mass = 5.9722e24                                          # kg
earth_angular_velocity = 2*np.pi/60/60/24                       # radians/sec earth rotation TODO: Get more accurate value

coordinates = [0.0, earth_radius+1000.0, 0.0]                   # m - note: coordinate convention is x,y,z global, with +z in the direction of global north
surface_rotation_rate = np.cross([0,0,earth_angular_velocity],[coordinates[0],coordinates[1],0])
velocity = [surface_rotation_rate[0], 0.0, 0.0]                 # m/s
acceleration = [0.0,0.0,0.0]                                    # m/s^2
vehicle_body_frame_quaternion = [1/np.sqrt(2),0,0,1/np.sqrt(2)] # (w, x, y, z) quaternion convention
vehicle_body_frame_quaternion = [0.6531,0.27105,0.27105,0.6531] # (w, x, y, z) quaternion convention
vehicle_body_frame_rate = [0.0,0.0,0.0]                         # rad/s Tait-Bryan convention
roll_rate = 0.0
pitch_rate = 0.0
yaw_rate = 0.0

# Instantiate Debugging Data Arrays
t = []

thrust_force_x_log = []
thrust_force_y_log = []
thrust_force_z_log = []

drag_force_x_log = []
drag_force_y_log = []
drag_force_z_log = []

gravity_force_x_log = []
gravity_force_y_log = []
gravity_force_z_log = []

coordinate_x_log = []
coordinate_y_log = []
coordinate_z_log = []

velocity_x_log = []
velocity_y_log = []
velocity_z_log = []

acceleration_x_log = []
acceleration_y_log = []
acceleration_z_log = []

stage_flag = 0

# Run Simulation
for i in np.linspace(0,simulation_time,int(simulation_time/dt)):
    distance_to_earth_center = np.linalg.norm(coordinates)
    if distance_to_earth_center < earth_radius:
        print('Rocket has landed rapidly (crashed) at time:',i)
        break

    # Get atmospheric parameters
    pressure,temperature,density = atmospheric_model.get_params(distance_to_earth_center)

    # Update vehicle parameters
    vehicle_mass = vehicle_model.getMass()
    drag_cross_section = vehicle_model.cross_section_area

    # Calculate Gravity
    gravitational_force = G*earth_mass*vehicle_mass/(distance_to_earth_center**2)
    gravity_vector = np.negative(coordinates)/distance_to_earth_center
    gravity_force_x = gravity_vector[0]*gravitational_force
    gravity_force_y = gravity_vector[1]*gravitational_force
    gravity_force_z = gravity_vector[2]*gravitational_force

    # Calculate Drag
    free_stream_velocity = np.cross([0,0,earth_angular_velocity],[coordinates[0],coordinates[1],0])
    air_speed_vector = velocity - free_stream_velocity
    air_speed = np.linalg.norm(air_speed_vector)
    drag_coefficient = vehicle_model.getDragCoefficient(air_speed,density,pressure)
    drag_force = 0.5*density*drag_cross_section*drag_coefficient*air_speed**2
    velocity_normalization = np.linalg.norm(velocity)
    drag_force_x = -drag_force*velocity[0]/velocity_normalization
    drag_force_y = -drag_force*velocity[1]/velocity_normalization
    drag_force_z = -drag_force*velocity[2]/velocity_normalization

    # Calculate Thrust
    thrust_force, thrust_vector, engine_distance = vehicle_model.getThrustVector(pressure,vehicle_body_frame_quaternion,i,dt) # TODO: future-proof this for control law
    thrust_force_x = thrust_force*thrust_vector[0]
    thrust_force_y = thrust_force*thrust_vector[1]
    thrust_force_z = thrust_force*thrust_vector[2]

    # Sum forces
    forces_x = gravity_force_x + drag_force_x + thrust_force_x
    forces_y = gravity_force_y + drag_force_y + thrust_force_y
    forces_z = gravity_force_z + drag_force_z + thrust_force_z

    # Calculate Velocities
    acceleration[0] = forces_x/vehicle_mass
    acceleration[1] = forces_y/vehicle_mass
    acceleration[2] = forces_z/vehicle_mass

    velocity[0] += acceleration[0]*dt
    velocity[1] += acceleration[1]*dt
    velocity[2] += acceleration[2]*dt
    
    # Convert to body coordinates - use quaternions
    global_to_body_quaternion = [vehicle_body_frame_quaternion[0],-vehicle_body_frame_quaternion[1],-vehicle_body_frame_quaternion[2],-vehicle_body_frame_quaternion[3]]
    rotated_thrust_force = general_functions.vector_quaternion_rotation([thrust_force_x,thrust_force_y,thrust_force_z],global_to_body_quaternion)
    rotated_drag_force = general_functions.vector_quaternion_rotation([drag_force_x,drag_force_y,drag_force_z],global_to_body_quaternion)
    rotated_gravity_force = general_functions.vector_quaternion_rotation([gravity_force_x,gravity_force_y,gravity_force_z],global_to_body_quaternion)
    # print(rotated_thrust_force)
    # Calculate Moments
    '''Moment distances are considered as positive values, measured as distance from payload (satellite)'''
    center_of_gravity = vehicle_model.getGravityCenter()
    center_of_pressure = vehicle_model.getCenterofPressure()
    center_of_pressure = center_of_gravity # TODO: Remove this, assumes drag acts with no moment on the rocket. Totally wrong but useful for validation of pure translation math
    Ixx, Iyy, Izz = vehicle_model.getMOIs()

    engine_cg_distance = center_of_gravity - engine_distance
    center_of_pressure_cg_distance = center_of_gravity - center_of_pressure

    '''Uses T = r x F to produce moment vector of form [r - roll,p - pitch,q - yaw]'''
    moments_vector = np.cross([-engine_cg_distance,0,0],rotated_thrust_force) + np.cross([-center_of_pressure_cg_distance,0,0],rotated_drag_force)

    # Calculate angular rates
    roll_rate += moments_vector[0]/Izz*dt
    pitch_rate += moments_vector[1]/Iyy*dt
    yaw_rate += moments_vector[2]/Ixx*dt

    # Update Step
    coordinates[0] += velocity[0]*dt 
    coordinates[1] += velocity[1]*dt 
    coordinates[2] += velocity[2]*dt
    
    vehicle_body_frame_quaternion = general_functions.quaternion_rotate_by_rate(vehicle_body_frame_quaternion,roll_rate,pitch_rate,yaw_rate,dt)

    if stage_flag != vehicle_model.stage_flag:
        if stage_flag == 1:
            stage_one_end_coordinates = list(coordinates)
        else:
            stage_two_end_coordinates = list(coordinates)
        stage_flag = vehicle_model.stage_flag
        print('Stage ', stage_flag, ' commenced at time T+', round(i,3))

    # Datalogging
    t.append(i)

    thrust_force_x_log.append(thrust_force_x)
    thrust_force_y_log.append(thrust_force_y)
    thrust_force_z_log.append(thrust_force_z)

    drag_force_x_log.append(drag_force_x)
    drag_force_y_log.append(drag_force_y)
    drag_force_z_log.append(drag_force_z)

    gravity_force_x_log.append(gravity_force_x)
    gravity_force_y_log.append(gravity_force_y)
    gravity_force_z_log.append(gravity_force_z)

    coordinate_x_log.append(coordinates[0])
    coordinate_y_log.append(coordinates[1])
    coordinate_z_log.append(coordinates[2])

    velocity_x_log.append(velocity[0])
    velocity_y_log.append(velocity[1])
    velocity_z_log.append(velocity[2])

    acceleration_x_log.append(acceleration[0])
    acceleration_y_log.append(acceleration[1])
    acceleration_z_log.append(acceleration[2])

    # print("Thrust:",thrust_force,"\tGravity:",gravitational_force,"\tDrag:",drag_force,"\tMass:",vehicle_mass)
    # print(free_stream_velocity, distance_to_earth_center-earth_radius)



plt.style.use('dark_background')
ax = plt.figure().subplot_mosaic(
    [["globe","globe","coordinate_x","velocity_x","acceleration_x"],
     ["globe","globe","coordinate_y","velocity_y","acceleration_y"],
     ["globe","globe","coordinate_z","velocity_z","acceleration_z"],
     ["globe","globe","thrust_x","drag_x","gravity_x"],
     ["globe","globe","thrust_y","drag_y","gravity_y"],
     ["globe","globe","thrust_z","drag_z","gravity_z"]],
     per_subplot_kw={"globe": {"projection": "3d"}}
)
# Configure 3-D Plot
frame_multiplier = 1.5
frame_limit = frame_multiplier*earth_radius
ax["globe"].set_xlim([-frame_limit,frame_limit])
ax["globe"].set_ylim([-frame_limit,frame_limit])
ax["globe"].set_zlim([-frame_limit,frame_limit])
ax["globe"].set_box_aspect([1,1,1])
ax["globe"].plot(coordinate_x_log,coordinate_y_log,coordinate_z_log)
# ax["globe"].plot(stage_one_end_coordinates[0],stage_one_end_coordinates[1],stage_one_end_coordinates[2],"r",marker = 'o',markersize=1)
# ax["globe"].plot(stage_two_end_coordinates[0],stage_two_end_coordinates[1],stage_two_end_coordinates[2],"r",marker = 'o',markersize=1)
ax["globe"].set_xlabel('X [m]')
ax["globe"].set_ylabel('Y [m]')
ax["globe"].set_zlabel('Z [m]')

# Generate Translucent Earth
u, v = np.mgrid[0:2*np.pi:64j, 0:np.pi:64j]
x_earth = earth_radius*np.cos(u)*np.sin(v)
y_earth = earth_radius*np.sin(u)*np.sin(v)
z_earth = earth_radius*np.cos(v)
ax["globe"].plot_surface(x_earth, y_earth, z_earth, color="b",alpha=0.25)

# Create Graphs
ax["coordinate_x"].plot(t,coordinate_x_log)
ax["coordinate_x"].set_title('Coordinates')
ax["coordinate_x"].set_xlabel('time [s]')
ax["coordinate_x"].set_ylabel('$X$ [m]')
ax["coordinate_y"].plot(t,coordinate_y_log)
ax["coordinate_y"].set_xlabel('time [s]')
ax["coordinate_y"].set_ylabel('$Y$ [m]')
ax["coordinate_z"].plot(t,coordinate_z_log)
ax["coordinate_z"].set_xlabel('time [s]')
ax["coordinate_z"].set_ylabel('$Z$ [m]')
ax["velocity_x"].plot(t,velocity_x_log)
ax["velocity_x"].set_title('Velocities')
ax["velocity_x"].set_xlabel('time [s]')
ax["velocity_x"].set_ylabel('$V_X$ [N]')
ax["velocity_y"].plot(t,velocity_y_log)
ax["velocity_y"].set_xlabel('time [s]')
ax["velocity_y"].set_ylabel('$V_Y$ [N]')
ax["velocity_z"].plot(t,velocity_z_log)
ax["velocity_z"].set_xlabel('time [s]')
ax["velocity_z"].set_ylabel('$V_Z$ [N]')
ax["acceleration_x"].plot(t,acceleration_x_log)
ax["acceleration_x"].set_title('Accelerations')
ax["acceleration_x"].set_xlabel('time [s]')
ax["acceleration_x"].set_ylabel('$A_X$ [N]')
ax["acceleration_y"].plot(t,acceleration_y_log)
ax["acceleration_y"].set_xlabel('time [s]')
ax["acceleration_y"].set_ylabel('$A_Y$ [N]')
ax["acceleration_z"].plot(t,acceleration_z_log)
ax["acceleration_z"].set_xlabel('time [s]')
ax["acceleration_z"].set_ylabel('$A_Z$ [N]')
ax["thrust_x"].plot(t,thrust_force_x_log)
ax["thrust_x"].set_title('Thrust Forces')
ax["thrust_x"].set_xlabel('time [s]')
ax["thrust_x"].set_ylabel('$F_{T,X}$ [N]')
ax["thrust_y"].plot(t,thrust_force_y_log)
ax["thrust_y"].set_xlabel('time [s]')
ax["thrust_y"].set_ylabel('$F_{T,Y}$ [N]')
ax["thrust_z"].plot(t,thrust_force_z_log)
ax["thrust_z"].set_xlabel('time [s]')
ax["thrust_z"].set_ylabel('$F_{T,Z}$ [N]')
ax["drag_x"].plot(t,drag_force_x_log)
ax["drag_x"].set_title('Drag Forces')
ax["drag_x"].set_xlabel('time [s]')
ax["drag_x"].set_ylabel('$F_{D,X}$ [N]')
ax["drag_y"].plot(t,drag_force_y_log)
ax["drag_y"].set_xlabel('time [s]')
ax["drag_y"].set_ylabel('$F_{D,Y}$ [N]')
ax["drag_z"].plot(t,drag_force_z_log)
ax["drag_z"].set_xlabel('time [s]')
ax["drag_z"].set_ylabel('$F_{D,Z}$ [N]')
ax["gravity_x"].plot(t,gravity_force_x_log)
ax["gravity_x"].set_title('Gravity Forces')
ax["gravity_x"].set_xlabel('time [s]')
ax["gravity_x"].set_ylabel('$F_{G,X}$ [N]')
ax["gravity_y"].plot(t,gravity_force_y_log)
ax["gravity_y"].set_xlabel('time [s]')
ax["gravity_y"].set_ylabel('$F_{G,Y}$ [N]')
ax["gravity_z"].plot(t,gravity_force_z_log)
ax["gravity_z"].set_xlabel('time [s]')
ax["gravity_z"].set_ylabel('$F_{G,Z}$ [N]')

figManager = plt.get_current_fig_manager()
# figManager.window.showMaximized()
plt.subplots_adjust(wspace=0.5,hspace=1.25)
plt.show()