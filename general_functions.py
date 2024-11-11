import numpy as np
# https://stengel.mycpanel.princeton.edu/Quaternions.pdf

def angular_velocities_quaternion(roll_rate, pitch_rate, yaw_rate, q):
    omega_matrix = [
        [         0,    yaw_rate, -pitch_rate,  roll_rate],
        [ -yaw_rate,           0,   roll_rate, pitch_rate],
        [pitch_rate,  -roll_rate,           0,          0],
        [-roll_rate, -pitch_rate,   -yaw_rate,          0]
    ]
    q_rate = 1/2 * np.matmul(omega_matrix,q)
    return q_rate

def vector_quaternion_rotation(v,q):
    v_false_quaternion =  v + [0]
    q = q / np.linalg.norm(q)
    q_row = [
        [q[0]**2+q[1]**2-q[2]**2-q[3]**2,         2*q[1]*q[2]+2*q[0]*q[3],         2*q[1]*q[3]-2*q[0]*q[2],                               0],
        [        2*q[1]*q[3]-2*q[0]*q[3], q[0]**2-q[1]**2+q[2]**2-q[3]**2,         2*q[2]*q[3]+2*q[0]*q[1],                               0],
        [        2*q[1]*q[3]+2*q[0]*q[2],         2*q[2]*q[3]-2*q[0]*q[1], q[0]**2-q[1]**2-q[2]**2+q[3]**2,                               0],
        [                              0,                               0,                               0, q[0]**2+q[1]**2+q[2]**2+q[3]**2],
    ]
    v_rotated = np.matmul(v_false_quaternion,q_row)
    return v_rotated[:-1]

def quaternion_rotate_by_rate(quaternion,wx,wy,wz,dt):
    rotation_update_matrix = [
        [2/dt,  -wx,  -wy,  -wz],
        [  wx, 2/dt,   wz,  -wy],
        [  wy,  -wz, 2/dt,   wx],
        [  wz,   wy,  -wx, 2/dt]
        ]
    
    rotated_quaternion = dt/2 * np.matmul(rotation_update_matrix, quaternion)
    return rotated_quaternion/np.linalg.norm(rotated_quaternion)

def quaternion_to_euler_321(q):
    roll = np.arctan2((2*q[0]*q[1]+q[2]*q[3]),(1-2*(q[1]**2 + q[2]**2)))
    pitch = -np.pi/2 + 2*np.arctan2(np.sqrt(1+2*(q[0]*q[2] - q[1]*q[3])),np.sqrt(1-2*(q[0]*q[2] - q[1]*q[3])))
    yaw = np.arctan2((2*q[0]*q[3]+q[2]*q[1]),(1-2*(q[2]**2 + q[3]**2)))
    return roll, pitch, yaw

def angle_between_vectors(v1,v2):
    dot_product = np.dot(v1,v2)
    n1 = np.linalg.norm(v1)
    n2 = np.linalg.norm(v2)
    angle = np.arccos(dot_product/(n1*n2))
    return angle