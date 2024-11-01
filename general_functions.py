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