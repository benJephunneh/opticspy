import numpy as np
from numpy import cos as cos
from numpy import sin as sin
from numpy import arctan2 as arctan2
from numpy import sqrt as sqrt
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from . import tools as tools


class Coefficient(object):
	"""
	Return a set of Seidel wavsfront aberrations Coefficient
	"""
	coefficients = []
	seidellist_=["Ap Piston",
					"At Tilt",
					"Ad Defocus",
		 			"As Astigmatism",
		 			"Ac Coma",
		 			"As Spherical"]
	def __init__(self,Ap=0,Bp=0,At=0,Bt=0,Ad=0,Bd=0,Aa=0,Ba=0,Ac=0,Bc=0,As=0,Bs=0):
		if type(Ap) == list:
			self.coefficients = Ap
		else:
			self.coefficients = [[Ap,Bp],[At,Bt],[Ad,Bd],[Aa,Ba],[Ac,Bc],[As,Bs]]
	def outputcoefficient(self):
		return self.coefficients

		
	def seidelsurface(self, label = True, zlim=[], matrix = False):
		r1 = np.linspace(0, 1, 100)
		u1 = np.linspace(0, 2*np.pi, 100)
		[u,r] = np.meshgrid(u1,r1)
		X = r*cos(u)
		Y = r*sin(u)
		W = seidelpolar(self.coefficients,r,u)
		fig = plt.figure(figsize=(12, 8), dpi=80)
		ax = fig.gca(projection='3d')
		surf = ax.plot_surface(X, Y, W, rstride=1, cstride=1, cmap=cm.RdYlGn,
	        linewidth=0, antialiased=False, alpha = 0.6)
		fig.colorbar(surf, shrink=1, aspect=30)
		plt.show()



	def twyman_green(self, lambda_1 = 632, PR = 1):
		lambda_1 = lambda_1*(10**-9)
		A = self.coefficients
		r = np.linspace(-PR, PR, 400)
		x, y = np.meshgrid(r,r) 
		OPD = seidelcartesian(A,x,y)*2/PR
		ph = 2 * np.pi * OPD
		I1 = 1
		I2 = 1
		Ixy = I1 + I2 + 2 * np.sqrt(I1*I2) * np.cos(ph)
		tools.makecircle(Ixy, r, PR) 
		fig = plt.figure(figsize=(9, 6), dpi=80)
		plt.imshow(-Ixy, extent=[-PR,PR,-PR,PR])
		plt.set_cmap('Greys')
		plt.show()


def seidelpolar(coefficient,r,u):
	W = coefficient
	h = 1
	Ap = W[0][0] * h**2
	At = W[1][0] * h*r*cos(u-W[1][1]*np.pi/180)
	Ad = W[2][0] * r**2
	Aa = W[3][0] * h**2*r**2*cos(u-W[3][1]*np.pi/180)
	Ac = W[4][0] * h*r*cos(u-W[4][1]*np.pi/180)
	As = W[5][0] * r**4

	W = Ap+At+Ad+Aa+Ac+As
	return W

def seidelcartesian(coefficient,x,y):
	W = coefficient
	h = 1
	u = arctan2(y,x)
	r = sqrt(x**2+y**2)
	W = seidelpolar(coefficient,r,u)
	return W
