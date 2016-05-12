import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

r, rho, y, z = np.loadtxt('out.dat', delimiter='\t', unpack=True)
plt.plot(x,z, label='Loaded from file!')

plt.xlabel('Radius r')
plt.ylabel('Density rho')
#plt.title('Interesting Graph\nCheck it out').
plt.legend()
plt.show()
plt.savefig('/out/rho.png')


mpl.rcParams['legend.fontsize'] = 10

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x, y, z, label='parametric curve')
ax.legend()

plt.show()
