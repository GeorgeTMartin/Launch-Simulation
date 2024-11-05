import numpy as np
# https://stengel.mycpanel.princeton.edu/Quaternions.pdf

def coordinateTransform(x,twist):
    gamma = twist[0]
    beta = twist[1]
    alpha = twist[2]
    rotation_matrix = [[np.cos(alpha)*np.cos(beta), np.cos(alpha)*np.sin(beta)*np.sin(gamma)-np.sin(alpha)*np.cos(gamma), np.cos(alpha)*np.sin(beta)*np.cos(gamma)+np.sin(alpha)*np.sin(gamma)],
                       [np.sin(alpha)*np.cos(beta), np.sin(alpha)*np.sin(beta)*np.sin(gamma)+np.cos(alpha)*np.cos(gamma), np.sin(alpha)*np.sin(beta)*np.cos(gamma)-np.cos(alpha)*np.sin(gamma)],
                       [-np.sin(beta), np.cos(beta)*np.sin(gamma), np.cos(beta)*np.cos(gamma)]]
    rotated_vector = np.matmul(rotation_matrix, np.array(x))
    return rotated_vector

def inverseCoordinateTransform(x,twist):
    gamma = twist[0]
    beta = twist[1]
    alpha = twist[2]
    rotation_matrix = [[np.cos(alpha)*np.cos(beta), np.cos(alpha)*np.sin(beta)*np.sin(gamma)-np.sin(alpha)*np.cos(gamma), np.cos(alpha)*np.sin(beta)*np.cos(gamma)+np.sin(alpha)*np.sin(gamma)],
                       [np.sin(alpha)*np.cos(beta), np.sin(alpha)*np.sin(beta)*np.sin(gamma)+np.cos(alpha)*np.cos(gamma), np.sin(alpha)*np.sin(beta)*np.cos(gamma)-np.cos(alpha)*np.sin(gamma)],
                       [-np.sin(beta), np.cos(beta)*np.sin(gamma), np.cos(beta)*np.cos(gamma)]]
    inverse_matrix = np.linalg.inv(rotation_matrix)    
    globalized_vector = np.matmul(inverse_matrix, x)

    return globalized_vector

def angleRates_BodytoEuler(x, twist):
    # https://liqul.github.io/blog/assets/rotation.pdf
    theta = twist[0]
    phi = twist[1]
    # body_rate_to_euler_rate_matrix = [[1, np.sin(phi)*np.tan(theta), np.cos(phi)*np.tan(theta)],
    #                                   [0, np.cos(phi), -np.sin(phi)],
    #                                   [0, np.sin(phi)/np.cos(theta), np.cos(phi)/np.cos(theta)]]
    # euler_angle_rates = np.matmul(body_rate_to_euler_rate_matrix,x)

    # test_matrix = [[1,0,-np.sin(theta)],
    #                [0, np.cos(phi), np.sin(phi)*np.cos(theta)],
    #                [0, -np.sin(phi), np.cos(phi)*np.cos(theta)]]
    # test_inversion_matrix = np.linalg.inv(test_matrix)
    # euler_angle_rates = np.matmul(test_inversion_matrix, x)

    body_rate_to_euler_rate_matrix = [[0,np.sin(phi),np.cos(phi)],
                               [0, np.cos(phi)*np.cos(theta), -np.sin(phi)*np.cos(theta)],
                               [np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi)*np.sin(theta)]]
    euler_angle_rates = np.matmul(body_rate_to_euler_rate_matrix,x)/np.cos(theta)

    return euler_angle_rates