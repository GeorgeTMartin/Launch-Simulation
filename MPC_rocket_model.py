'''
Model References:
Longuski, J., Guzmán, J. J., & Prussing, J. E. (2014). Optimal control with aerospace applications. Springer. 

Section 7.3 - Time Optimal Launching of a Satellite
Assumes flat Earth model and simplified form to determine launch profile for reference in control law.
Desire altitude varying velocity vector (e.g. magnitude and direction).
'''
import matplotlib.pyplot as plt
from casadi import *
G = 6.67430*10**-11                                             # Earth Parameters - Point Mass Model
earth_radius = 6378137                                          # meters
earth_mass = 5.9722e24                                          # kg

orbit_altitude = 1.0e9                                          # meters
orbit_speed = np.sqrt(G*earth_mass/orbit_altitude)              # m/s 


def rocket_model():
    # define structs
    constraint = types.SimpleNamespace()
    model = types.SimpleNamespace()

    model_name = "FalconIX"
    ## Race car parameters
    m_s1_solid = 25600.0    # kg
    m_s1_liquid = 415300.0  # kg
    m_s2_solid = 2900.0     # kg
    m_s2_liquid = 92670.0   # kg
    m_payload = 12000.0     # kg
    m_dot = -2100.0         # kg/s
    exit_velocity = 3000.0  # m/s

    s1m = m_s1_solid+m_s1_liquid+m_s2_solid+m_s2_liquid+m_payload
    s2m = m_s2_solid+m_s2_liquid+m_payload

    t0 = 0                  # launch [s]
    t1 = 196.036            # s
    t2 = 591.107            # s

    # Piecewise-linear - derived from https://github.com/casadi/casadi/issues/2489
    t = SX.sym("t")
    stage_one_mass_function = s1m + m_dot*t             # Thrust of 9 engines
    stage_two_mass_function = s2m + (m_dot/9.0)*(t-t1)  # Thrust of 1 engine
    # stage_three_mass_function = SX(m_payload)           # No Thrust - Mass just SX(m_payload)

    times = SX.sym('times',6)
    times[0] = t0
    times[1] = t1-.0001       # Small margin to model rapid loss of stage solid weight
    times[2] = t1
    times[3] = t2-.0001
    times[4] = t2
    times[5] = t2+1000        # Making arbitrary range of const mass
    rocket_thrust_equation = SX.sym('rocket_mass',6)
    rocket_thrust_equation[0] = -m_dot*exit_velocity/(9.8*stage_one_mass_function)
    rocket_thrust_equation[1] = -m_dot*exit_velocity/(9.8*stage_one_mass_function)
    rocket_thrust_equation[2] = -m_dot*exit_velocity/(9.8*stage_two_mass_function)
    rocket_thrust_equation[3] = -m_dot*exit_velocity/(9.8*stage_two_mass_function)
    rocket_thrust_equation[4] = SX(0)
    rocket_thrust_equation[5] = SX(0)

    # rocket_thrust_function = Function('func',[t],[pw_lin(t,times,rocket_thrust_equation)])


    ## CasADi Model
    # set up states & controls
    x = SX.sym("x")
    y = SX.sym("y")
    x_dot = SX.sym("x_dot")
    y_dot = SX.sym("y_dot")
    theta = SX.sym("theta")
    theta_dot = SX.sym("theta_dot")

    states = vertcat(x, y, x_dot, y_dot, theta, t)

    # controls
    u = vertcat(theta_dot)

    # State Derivatives
    x_d_dot = SX.sym("x_d_dot")
    y_d_dot = SX.sym("y_d_dot")
    t_dot = SX.sym("t_dot")
    states_dot = vertcat(x_dot, y_dot, x_d_dot, y_d_dot, theta_dot, t_dot)

    # algebraic variables
    z = vertcat([])

    # parameters
    p = vertcat([])

    # dynamics
    Ispm = exit_velocity/9.8
    # thrust_acceleration = Function('thrust',[t],[])
    g = G*earth_mass/(y+earth_radius)**2
    # x_accel_function = Function('rocket_thrust_x',[t,theta],[rocket_thrust_equation*np.cos(theta)]) 
    # y_accel_function = Function('rocket_thrust_y',[t,theta,y],[rocket_thrust_equation*np.sin(theta)-g])  

    # Explicit Dynamics [x_dot = f(x,u)]
    f_expl = vertcat(
        x_dot,                                                  # X_dot
        y_dot,                                                  # Y_dot
        pw_lin(t,times,rocket_thrust_equation)*np.cos(theta),   # x_d_dot
        pw_lin(t,times,rocket_thrust_equation)*np.sin(theta)-g, # y_d_dot
        theta_dot,                                              # theta_dot
        1                                                       # t_dot [s/s]
    )


    # state bounds
    model.theta_min = 0.0  # minimum pitch angle [rad]
    model.theta_max = np.pi/2  # maximum pitch angle [rad]
    model.theta_dot_min = -0.1 # rad/s
    model.theta_dot_max = 0.1  # rad/s

    # Define initial conditions
    model.x0 = np.array([0, 0, 0, 0, np.pi/2, 0])

    # define constraints struct
    constraint.expr = vertcat(theta)

    # Define model struct
    params = types.SimpleNamespace()
    params.m  = rocket_thrust_equation
    params.x_dot = x_dot
    model.f_impl_expr = states_dot - f_expl
    model.f_expl_expr = f_expl
    model.x = states
    model.xdot = states_dot
    model.u = u
    model.z = z
    model.p = p
    model.name = model_name
    model.params = params
    return model, constraint



if __name__ == "__main__":
    rocket_model()