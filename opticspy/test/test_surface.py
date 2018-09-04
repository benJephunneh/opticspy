import numpy as np
from numpy import sqrt as sqrt
from numpy import cos as cos
from numpy import sin as sin
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib.ticker import LinearLocator as LinearLocator
from matplotlib.ticker import FormatStrFormatter as FormatStrFormatter
#generate test surface figure
def spherical_surf(l1):
	R = 1.02
	l1 = l1  #surface matrix length
	theta = np.linspace(0, 2*np.pi, l1)
	rho = np.linspace(0, 1, l1)
	[u,r] = np.meshgrid(theta,rho)
	X = r*cos(u)
	Y = r*sin(u)
	Z = sqrt(R**2-r**2)-sqrt(R**2-1)
	v_1 = max(abs(Z.max()),abs(Z.min()))

	noise = (np.random.rand(len(Z),len(Z))*2-1)*0.05*v_1
	Z = Z+noise
	fig = plt.figure(figsize=(12, 8), dpi=80)
	ax = fig.gca(projection='3d')
	surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.RdYlGn,\
								linewidth=0, antialiased=False, alpha = 0.6)
	v = max(abs(Z.max()),abs(Z.min()))
	ax.set_zlim(-1, 2)
	ax.zaxis.set_major_locator(LinearLocator(10))
	ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
	cset = ax.contourf(X, Y, Z, zdir='z', offset=-1, cmap=cm.RdYlGn)
	fig.colorbar(surf, shrink=1, aspect=30)
	plt.title('Test Surface: Spherical surface with some noise',fontsize=16)
	plt.show()

	#Generate test surface matrix from a detector
	x = np.linspace(-1, 1, l1)
	y = np.linspace(-1, 1, l1)
	[X,Y] = np.meshgrid(x,y)
	Z = sqrt(R**2-(X**2+Y**2))-sqrt(R**2-1)+noise
	for i in range(len(Z)):
		for j in range(len(Z)):
			if x[i]**2+y[j]**2>1:
				Z[i][j]=0
	return Z