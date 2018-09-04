import numpy as np
from numpy import cos as cos
from numpy import sin as sin
import matplotlib.pyplot as plt
from matplotlib import cm as cm

class Coefficient(object):
	"""
	Return a set of Seidel wavsfront aberrations Coefficient
	"""
	coefficients = []
	seidellist_=["W200 piston",
					 "W111 tilt",
					 "W020 defocus",
					 "W040 spherical",
					 "W131 coma",
					 "W222 astigmatism",
					 "W220 field curvature",
					 "W311 distortion"]
	def __init__(self, h=0, W200=0,W111=0,W020=0,W040=0,W131=0,W222=0,W220=0,W311=0):
		if type(h) == list:
			self.coefficients = h
		else:
			self.coefficients = [h,W200,W111,W020,W040,W131,W222,W220,W311]
	def outputcoefficient(self):
		return self.coefficients
	def listcoefficient(self):
		"""
		------------------------------------------------
		listcoefficient():

		List the coefficient in Coefficient

		------------------------------------------------
		"""
		print("h="+str(self.coefficients[0]))
		for i in range(len(self.coefficients)-1):
			print(self.seidellist_[i][0:4]+"="+\
				str(self.coefficients[i+1])+self.seidellist_[i][4:])
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



def seidelpolar(coefficient,r,u):
	h = coefficient[0]
	W = coefficient
	W200 = W[1] * h**2
	W111 = W[2] * h*r*cos(u)
	W020 = W[3] * r**2
	W040 = W[4] * r**4
	W131 = W[5] * h*r**3*cos(u)
	W222 = W[6] * h**2*r**2*cos(u)
	W220 = W[7] * h**2*r**2
	W311 = W[8] * h**3*r*cos(u)
	W = W200+W111+W020+W040+W131+W222+W220+W311

	return W

