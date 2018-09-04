import numpy as np
from numpy import sqrt as sqrt
from numpy import cos as cos
from numpy import sin as sin
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib.ticker import LinearLocator as LinearLocator
from matplotlib.ticker import FormatStrFormatter as FormatStrFormatter
from numpy.fft import fftshift as fftshift
from numpy.fft import ifftshift as ifftshift
from numpy.fft import fft2 as fft2

def apershow(obj):
	obj = -abs(obj)
	plt.imshow(obj)
	plt.set_cmap('Greys')
	plt.show()

l1 = 100
#Generate test surface matrix from a detector
x = np.linspace(-1, 1, l1)
y = np.linspace(-1, 1, l1)
[X,Y] = np.meshgrid(x,y)
r = sqrt(X**2+Y**2)
Z = sqrt(14)*(8*X**4-8*X**2*r**2+r**4)*(6*r**2-5)
for i in range(len(Z)):
	for j in range(len(Z)):
		if x[i]**2+y[j]**2>1:
			Z[i][j]=0

fig = plt.figure(1)
ax = fig.gca(projection='3d')
surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.RdYlGn,
    linewidth=0, antialiased=False, alpha = 0.6)

v = max(abs(Z.max()),abs(Z.min()))
ax.set_zlim(-v*5, v*5)
cset = ax.contourf(X, Y, Z, zdir='z', offset=-v*5, cmap=cm.RdYlGn)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
fig.colorbar(surf, shrink=1, aspect=30)
plt.show()

d = 400
A = np.zeros([d,d])
A[d/2-49:d/2+51,d/2-49:d/2+51] = Z
plt.imshow(A)
plt.show()

abbe = np.exp(1j*2*np.pi*A)
for i in range(len(abbe)):
	for j in range(len(abbe)):
		if abbe[i][j]==1:
			abbe[i][j]=0
fig = plt.figure(2)
AP = abs(fftshift(fft2(fftshift(abbe))))**2
AP = AP/AP.max()
plt.imshow(AP)
plt.show()