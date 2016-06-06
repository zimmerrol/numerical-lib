import matplotlib.pyplot as pyplot
import numpy as np

xmin = 0
xmax = 2
ymin = 0
xmax = 1.4
xres = 80
yres = 56

values = np.genfromtxt(fname = "res/out.dat", delimiter = '\t')
X = values[:,0]
Y = values[:,1]
Z = values[:,2]

print(len(Z))
print (Z)

pyplot.pcolor(Z.reshape(xres, yres).transpose(),)
pyplot.gca().set_aspect('equal', adjustable='box')
pyplot.draw()
pyplot.show()
