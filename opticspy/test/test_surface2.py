import numpy as np
from numpy import sqrt as sqrt
from numpy import cos as cos
from numpy import sin as sin
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib.ticker import LinearLocator as LinearLocator
from matplotlib.ticker import FormatStrFormatter as FormatStrFormatter
#generate test surface figure
def makecircle(a, r, PR):
	max = a.max()
	size = np.sqrt(a.size)
	for i in range(int(size)):
		for j in range(int(size)):
			if np.sqrt(r[i]**2+r[j]**2) > PR:
				a[i,j] = max

def testsurface2():

	lambda_1 = 632*(10**-9)
	PR = 1
	r = np.linspace(-PR, PR, 200)
	x, y = np.meshgrid(r,r) 
	r1 = np.sqrt(x**2 + y**2)
	Z4 = 1
	Z5 = 0.6
	ZX  =  Z4  * np.sqrt(3)*(2*r1**2-1) + Z5*2*np.sqrt(6)*x*y
	OPD = 	ZX*2/PR
	ph = 2 * np.pi * OPD
	Ia = 1
	Ib = 1
	Ixy = Ia + Ib + 2 * np.sqrt(Ia*Ib) * np.cos(ph)
	makecircle(Ixy, r, PR)
	fig = plt.figure(figsize=(9, 6), dpi=80)
	plt.imshow(-Ixy, extent=[-PR,PR,-PR,PR])
	plt.set_cmap('Greys')
	plt.show()

	I1 = Ia + Ib + 2 * np.sqrt(Ia*Ib) * np.cos(ph)
	I2 = Ia + Ib + 2 * np.sqrt(Ia*Ib) * np.cos(ph+45.0/180*np.pi)
	I3 = Ia + Ib + 2 * np.sqrt(Ia*Ib) * np.cos(ph+90.0/180*np.pi)
	I4 = Ia + Ib + 2 * np.sqrt(Ia*Ib) * np.cos(ph+135.0/180*np.pi)

	Ilist = [I1,I2,I3,I4]

	for i in range(4):
		makecircle(Ilist[i], r, PR)
		fig = plt.figure(figsize=(9, 6), dpi=80)
		plt.imshow(-Ilist[i], extent=[-PR,PR,-PR,PR])
		plt.set_cmap('Greys')
		plt.show()

	ph1 = np.arctan((I4-I2)/(I1-I3))

	Ixy1 = Ia + Ib + 2 * np.sqrt(Ia*Ib) * np.cos(ph1)
	fig = plt.figure(figsize=(9, 6), dpi=80)
	plt.imshow(-Ixy, extent=[-PR,PR,-PR,PR])
	plt.set_cmap('Greys')
	plt.show()

	OPD = ph*PR/2
	Z = OPD
	fig = plt.figure(figsize=(6, 6), dpi=80)
	#ax = fig.gca(projection='3d')
	#surf = ax.plot_surface(x, y, Z, rstride=1, cstride=1, cmap=cm.RdYlGn,linewidth=0, antialiased=False, alpha = 0.6)
	im = plt.pcolormesh(x, y, Z, cmap=cm.RdYlGn)
	plt.colorbar()
	plt.show()

	for i in range(len(Z)):
		for j in range(len(Z)):
			if r[i]**2+r[j]**2>1:
				Z[i][j]=0
	fig = plt.figure(figsize=(6, 6), dpi=80)
	im = plt.pcolormesh(x, y, Z, cmap=cm.RdYlGn)
	plt.colorbar()
	plt.show()

	return Z