import numpy as np
from numpy.linalg import inv
theta_max=16
points = 10
covariance_inverse = np.zeros((points,points))
covariance = np.loadtxt('covariance_matrix')[3:theta_max-3,3:theta_max-3]
covariance_inverse = inv(covariance)
np.savetxt('inv_matrix10',covariance_inverse,fmt='%.8f')