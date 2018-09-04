from __future__ import division as division
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm as cm

#from . import zernike as zernike
from . import interferometer_zenike as interferometer

from . import tools as tools

def hartmann(coefficients, r, R):
	"""
	Generate Hartmann spotdiagram
	use circle hartmann plate
	coefficients: zernike coefficients
	r: distance from the pupil of the wavefront to detect plate
	R: radius of mirror under test
	"""
	coefficients = coefficients.coefficients
	x_list = []
	y_list = []
	Ax_list = []
	Ay_list = []
	s = 40
	x = y = np.linspace(-R, R, s)
	M = []
	for i in x:
		tmp = []
		for j in y:
			if i**2 + j**2 < R**2:
				x_list.append(i)
				y_list.append(j)
				W0 = interferometer.zernikecartesian(coefficients,i,j)
				Wx = interferometer.zernikecartesian(coefficients,1.01*i,j)
				Wy = interferometer.zernikecartesian(coefficients,i,1.01*j)
				TAx = -(Wx-W0)/(0.01*i)*r
				TAy = -(Wy-W0)/(0.01*j)*r
				Ax_list.append(TAx)
				Ay_list.append(TAy)
				tmp.append([1,[i,j],[TAx,TAy]])
			else:
				tmp.append([0,[i,j],[0,0]])
		M.append(tmp)
	fig = plt.figure(1,figsize=(6, 6))
	ax = fig.gca()
	ax.set_axis_bgcolor('black')
	plt.title('Hartmann Spotdiagram',fontsize=18)
	plt.plot(Ax_list,Ay_list,'wo')
	plt.show()

	return M,r


def hartmann_rebuild(M,r):
	s = len(M)
	w = np.zeros([s,s])
	d = 2
	for n in range(s):
		label = 0
		for m in range(s):
			if M[n][m][0] == 0:
				pass
			elif (M[n][m][0] != 0 and label == 0):
				w[n,m] = 0
				label = 1
			elif (M[n][m][0] != 0 and label == 1):
				w[n,m] = w[n][m-1] + d/2/r*(M[n][m-1][2][0] + M[n][m][2][0])
			else:
				print('wrong')
	fig = plt.figure(2,figsize=(6, 6))
	plt.imshow(w)
	plt.show()
	# x = np.linspace(-1,1,s)
	# [X,Y] = np.meshgrid(x,x)
	# fig = plt.figure(figsize=(8, 8), dpi=80)
	# ax = fig.gca(projection='3d')
	# surf = ax.plot_surface(w, rstride=1, cstride=1, cmap=cm.RdYlGn,
	#         linewidth=0, antialiased=False, alpha = 0.6)
	# plt.show()

	return w


#Depth first search algorithm, use to find wavefrontase map(where)
def DFS(M,wavefront1,m,n,s):
	stack = []
	stack.append([m,n])
	M[m,n] = 2
	wavefront = np.zeros([s,s])

	while (len(stack) != 0):
		[m,n] = stack[-1]
		if m + 1 < s and n < s and M[m+1,n] == 1 and M[m+1,n] != 0 and M[m+1,n] != 2:
			m = m + 1
			M[m,n] = 2
			stack.append([m,n])

			wavefront[m,n] = wavefront[m-1,n] + v(wavefront1[m,n] - wavefront1[m-1,n])

		elif m - 1 > 0 and n < s and M[m-1,n] == 1 and M[m-1,n] != 0 and M[m-1,n] != 2:
			m = m - 1
			M[m,n] = 2
			stack.append([m,n])

			wavefront[m,n] = wavefront[m+1,n] + v(wavefront1[m,n] - wavefront1[m+1,n])

		elif m < s and n + 1 < s and M[m,n+1] == 1 and M[m,n+1] != 0 and M[m,n+1] != 2:
			n = n + 1
			M[m,n] = 2
			stack.append([m,n])

			wavefront[m,n] = wavefront[m,n-1] + v(wavefront1[m,n] - wavefront1[m,n-1])

		elif m < s and n - 1 > 0 and M[m,n-1] == 1 and M[m,n-1] != 0 and M[m,n-1] != 2:
			n = n - 1
			M[m,n] = 2
			stack.append([m,n])

			wavefront[m,n] = wavefront[m,n+1] + v(wavefront1[m,n] - wavefront1[m,n+1])

		else:
			stack.pop()
	return wavefront



