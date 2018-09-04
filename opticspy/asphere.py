from __future__ import division as division
import numpy as np
from numpy import sqrt as sqrt
from numpy import cos as cos
from numpy import sin as sin
import matplotlib.pyplot as plt
from matplotlib import cm as cm

class Coefficient(object):
	"""
	Return a set of Asphere Coefficient
	R,k,a2,a3,a4,a5,a6,a7,a8,a9,a10
	"""
	coefficients = []
	def __init__(self,R=0,k=0,a2=0,a3=0,a4=0,a5=0,a6=0,a7=0,a8=0,a9=0,a10=0):

		if type(R) == list:
			self.coefficients = R + [0]*(11-len(R))
		else:
			self.coefficients = [R,k,a2,a3,a4,a5,a6,a7,a8,a9,a10]
			
	def outputcoefficient(self):
		return self.coefficients

	def aspheresurface(self):
		"""
		Show the surface of an asphere.
		=============================================================
		Try: 
		A = opticspy.asphere.Coefficient(R=50,a2=0.18*10**(-8),a3 = 0.392629*10**(-13))

		"""
		R = self.coefficients[0]
		theta = np.linspace(0, 2*np.pi, 100)
		rho = np.linspace(0, R, 100)
		[u,r] = np.meshgrid(theta,rho)
		X = r*cos(u)
		Y = r*sin(u)
		Z = aspherepolar(self.coefficients,r)
		fig = plt.figure(figsize=(12, 8), dpi=80)
		ax = fig.gca(projection='3d')
		surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.RdYlGn,
	        linewidth=0, antialiased=False, alpha = 0.6)
		plt.show()
		return 0

	def aspherematrix(self):
		l = 100
		R = self.coefficients[0]
		x1 = np.linspace(-R, R, l)
		[X,Y] = np.meshgrid(x1,x1)
		r = sqrt(X**2+Y**2)
		Z = aspherepolar(self.coefficients,r)
		for i in range(l):
			for j in range(l):
				if x1[i]**2+x1[j]**2 > R**2:
					Z[i][j] = 0
		fig = plt.figure(figsize=(12, 8), dpi=80)
		ax = fig.gca(projection='3d')
		surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.RdYlGn,
	        linewidth=0, antialiased=False, alpha = 0.6)
		plt.show()
		return Z

	def asphereline(self):
		R,k,a2,a3,a4,a5,a6,a7,a8,a9,a10 = self.coefficients
		r = np.linspace(-R,R,100)
		C = 1/R
		Z = C*r**2*(1+sqrt(1-(1+k)*r**2*C**2)) + a2*r**4 + a3*r**6 + a4*r**8 + \
		+ a5*r**10 + a6*r**12 + a7*r**14 + a8*r**16 + a9*r**18 + a10*r**20
		Z = -Z
		fig = plt.figure(figsize=(12, 8), dpi=80)
		plt.plot(r,Z)
		plt.axis('equal')
		plt.show()

def aspherepolar(coefficient,r):
	R,k,a2,a3,a4,a5,a6,a7,a8,a9,a10 = coefficient
	C = 1/R
	Z = C*r**2*(1+sqrt(1-(1+k)*r**2*C**2)) + a2*r**4 + a3*r**6 + a4*r**8 + \
		+ a5*r**10 + a6*r**12 + a7*r**14 + a8*r**16 + a9*r**18 + a10*r**20
	return -Z












