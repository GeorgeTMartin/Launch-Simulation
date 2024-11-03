import numpy as np


def coordinateTransform(x,twist):
    theta = twist[0]
    phi = twist[1]
    psi = twist[2]
    rotation_matrix = [[np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)],
                        [-np.cos(psi)*np.sin(phi)+np.sin(psi)*np.sin(theta)*np.cos(phi), np.cos(psi)*np.cos(phi)+np.sin(psi)*np.sin(theta)*np.sin(phi),np.sin(psi)*np.cos(theta)],
                        [np.sin(psi)*np.sin(theta)+np.cos(psi)*np.sin(theta)*np.cos(phi), -np.sin(psi)*np.cos(theta)+np.cos(psi)*np.sin(theta)*np.sin(phi), np.cos(psi)*np.cos(theta)]]
    rotated_vector = np.matmul(rotation_matrix, x)

    return rotated_vector

def inverseCoordinateTransform(x,twist):
    theta = twist[0]
    phi = twist[1]
    psi = twist[2]
    rotation_matrix = [[np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)],
                        [-np.cos(psi)*np.sin(phi)+np.sin(psi)*np.sin(theta)*np.cos(phi), np.cos(psi)*np.cos(phi)+np.sin(psi)*np.sin(theta)*np.sin(phi),np.sin(psi)*np.cos(theta)],
                        [np.sin(psi)*np.sin(theta)+np.cos(psi)*np.sin(theta)*np.cos(phi), -np.sin(psi)*np.cos(theta)+np.cos(psi)*np.sin(theta)*np.sin(phi), np.cos(psi)*np.cos(theta)]]
    inverse_matrix = np.linalg.inv(rotation_matrix)    
    globalized_vector = np.matmul(inverse_matrix, x)

    return globalized_vector

def angleRates_BodytoEuler(x, twist):
    theta = twist[0]
    phi = twist[1]
    inverse_rotation_matrix = [[0,np.sin(phi),np.cos(phi)],
                               [0, np.cos(phi)*np.cos(theta), -np.sin(phi)*np.cos(theta)],
                               [np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi)*np.sin(theta)]]
    euler_angle_rates = np.matmul(inverse_rotation_matrix,x)/np.cos(theta)
    return euler_angle_rates