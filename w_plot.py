import numpy as np
import matplotlib.pyplot as plt

theta_max = 16

x1 = np.logspace(-2.4,-0.6,10)
x2 = np.logspace(-3,0,16)

y1 = np.loadtxt('chain_best')
y2 = np.loadtxt('w_mean')

covariance_matrix  = np.loadtxt('covariance_matrix')
covariance_old     = covariance_matrix.reshape((theta_max,theta_max))
covariance 	       = covariance_old
bar = np.zeros(theta_max)
for i in range(theta_max):
	bar[i] = np.sqrt(covariance[i][i])
plt.errorbar(x2[2:], y2[2:theta_max], yerr=bar[2:theta_max], capsize=2, marker='s', markersize=1, color='r', label='w')
plt.plot(x1,y1,'c',label="HOD best-fit")
plt.xlabel("$\\theta$(deg)",size="xx-large")
plt.ylabel('w',size="xx-large")
plt.xscale('log')
plt.yscale('log',nonposy='clip')
plt.legend()
plt.show()