import numpy as np
import matplotlib.pyplot as plt
import corner
from scipy.stats import norm

ndim = 2
data = np.loadtxt('foutput')
data = data.reshape(7,10,10000)
draw = data[:,:,1000:]
sample = (draw[0:ndim,:,:].reshape([ndim,-1])).T

#value1 = np.array([50,1.0])
#value2 = np.mean(sample, axis=0)
a,b=np.where(draw[6]==np.min(draw[6]))
value2 = draw[0:ndim,a,b]


figure = corner.corner(sample,bins=[101,101],levels=(norm.cdf([1])),
labels=["$logM_{min}$", "$logM_{1}$"],plot_contours=True,
show_titles=False,title_kwargs={"fontsize": 15},label_kwargs={"fontsize": 20},title_fmt=".6f")
axes = np.array(figure.axes).reshape((ndim, ndim))
# Loop over the diagonal
#for i in range(ndim):
#    ax = axes[i, i]
#    ax.axvline(value1[i], color="g")
#    ax.axvline(value2[i], color="r")

# Loop over the histograms
for yi in range(ndim):
    for xi in range(yi):
        ax = axes[yi, xi]
#        ax.axvline(value1[xi], color="g")
#        ax.axvline(value2[xi], color="r")
#        ax.axhline(value1[yi], color="g")
#        ax.axhline(value2[yi], color="r")
#        ax.plot(value1[xi], value1[yi], "sg")
        ax.plot(value2[xi], value2[yi], "xr",markersize=10)

#sigam_best = np.percentile(sample[:,0],50)
#sigam_less = np.percentile(sample[:,0],16)
#sigma_more = np.percentile(sample[:,0],84)
#bias_best = np.percentile(sample[:,1],50)
#bias_less = np.percentile(sample[:,1],16)
#bias_more = np.percentile(sample[:,1],84)
#print(sigam_best,sigam_less,sigma_more,bias_best,bias_less,bias_more)
print(value2)
plt.show()