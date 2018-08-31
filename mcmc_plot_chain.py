import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(1,10000,10000)

data = np.loadtxt('foutput')

data = data.reshape(7,10,10000)

plt.subplot(421)
plt.plot(x,data[0,5])
plt.title('$logM_{min}$')

plt.subplot(422)
plt.plot(x,data[1,5])
plt.title('$logM_{1}$')

plt.subplot(423)
plt.plot(x,data[2,5])
plt.title('$n_{g}$')

plt.subplot(424)
plt.plot(x,data[3,5])
plt.title('$M_{h}$')

plt.subplot(425)
plt.plot(x,data[4,5])
plt.title('$bias$')

plt.subplot(426)
plt.plot(x,data[5,5])
plt.title('$f_{sat}$')

plt.subplot(427)
plt.plot(x,data[6,5])
plt.title('${\chi}^{2}$')
plt.show()