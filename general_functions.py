import numpy as np
# https://stengel.mycpanel.princeton.edu/Quaternions.pdf

def angular_velocities_quaternion(q1, q2, dt):
    return (2 / dt) * np.array([
        q1[0]*q2[1] - q1[1]*q2[0] - q1[2]*q2[3] + q1[3]*q2[2],
        q1[0]*q2[2] + q1[1]*q2[3] - q1[2]*q2[0] - q1[3]*q2[1],
        q1[0]*q2[3] - q1[1]*q2[2] + q1[2]*q2[1] - q1[3]*q2[0]])

def vector_quaternion_rotation(v,q):
    v_quaternion = v + [0]
    q = q / np.linalg.norm(q)
    q_row = [
        [q[0]**2+q[1]**2-q[2]**2-q[3]**2,         2*q[1]*q[2]+2*q[0]*q[3],         2*q[1]*q[3]-2*q[0]*q[2],                               0],
        [        2*q[1]*q[3]-2*q[0]*q[3], q[0]**2-q[1]**2+q[2]**2-q[3]**2,         2*q[2]*q[3]+2*q[0]*q[1],                               0],
        [        2*q[1]*q[3]+2*q[0]*q[2],         2*q[2]*q[3]-2*q[0]*q[1], q[0]**2-q[1]**2-q[2]**2+q[3]**2,                               0],
        [                              0,                               0,         2*q[1]*q[3]-2*q[0]*q[2], q[0]**2+q[1]**2+q[2]**2+q[3]**2],
    ]
    v_rotated = np.matmul(v_quaternion,q_row)
    return v_rotated[0:3]