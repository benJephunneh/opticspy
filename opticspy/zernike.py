# coding=utf-8
from __future__ import division as division
import numpy as np
from numpy import cos as cos
from numpy import sin as sin
from numpy import sqrt as sqrt
from numpy import arctan2 as arctan2
import math
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib.ticker import LinearLocator as LinearLocator
from matplotlib.ticker import FormatStrFormatter as FormatStrFormatter
from numpy.fft import fftshift as fftshift
from numpy.fft import ifftshift as ifftshift
from numpy.fft import fft2 as fft2
from numpy.linalg import inv
from numpy.fft import ifft2 as ifft2
import operator as op

from . import interferometer_zenike as interferometer
from . import seidel2 as seidel2
from . import tools as tools
from . import hartmann as hartmann

class Coefficient(object):
	"""
	Return a set of Zernike Polynomials Coefficient
	"""
	coefficients = []
	zernikelist = [ "Z00 Piston or Bias",
						"Z11 x Tilt",
						"Z11 y Tilt",
						"Z20 Defocus",
						"Z22 Primary Astigmatism at 45",
						"Z22 Primary Astigmatism at 0",
						"Z31 Primary y Coma",
						"Z31 Primary x Coma",
						"Z33 y Trefoil",
						"Z33 x Trefoil",
						"Z40 Primary Spherical",
						"Z42 Secondary Astigmatism at 0",
						"Z42 Secondary Astigmatism at 45",
						"Z44 x Tetrafoil",
						"Z44 y Tetrafoil",
						"Z51 Secondary x Coma",
						"Z51 Secondary y Coma",
						"Z53 Secondary x Trefoil",
						"Z53 Secondary y Trefoil",
						"Z55 x Pentafoil",
						"Z55 y Pentafoil",
						"Z60 Secondary Spherical",
						"Z62 Tertiary Astigmatism at 45",
						"Z62 Tertiary Astigmatism at 0",
						"Z64 Secondary x Trefoil",
						"Z64 Secondary y Trefoil",
						"Z66 Hexafoil Y",
						"Z66 Hexafoil X",
						"Z71 Tertiary y Coma",
						"Z71 Tertiary x Coma",
						"Z73 Tertiary y Trefoil",
						"Z73 Tertiary x Trefoil",
						"Z75 Secondary Pentafoil Y",
						"Z75 Secondary Pentafoil X",
						"Z77 Heptafoil Y",
						"Z77 Heptafoil X",
						"Z80 Tertiary Spherical"]

	def __init__(self,
			Z1=0, Z2=0, Z3=0, Z4=0, Z5=0, Z6=0, Z7=0, \
			Z8=0, Z9=0, Z10=0, Z11=0, Z12=0, Z13=0, Z14=0, \
			Z15=0, Z16=0, Z17=0, Z18=0, Z19=0, Z20=0, Z21=0, \
			Z22=0, Z23=0, Z24=0, Z25=0, Z26=0, Z27=0, Z28=0, \
			Z29=0, Z30=0, Z31=0, Z32=0, Z33=0, Z34=0, Z35=0, Z36=0, Z37=0):
		if type(Z1) == list:
			self.coefficients = Z1 + [0]*(37-len(Z1))
		else:
			self.coefficients = [Z1, Z2, Z3, Z4, Z5, Z6, Z7,
					Z8, Z9, Z10, Z11, Z12, Z13, Z14, Z15, Z16, Z17,
					Z18, Z19, Z20, Z21, Z22, Z23, Z24, Z25, Z26,
					Z27, Z28, Z29, Z30, Z31, Z32, Z33, Z34, Z35, Z36, Z37]
	def outputcoefficient(self):
		return self.coefficients
	def listcoefficient(self):
		"""
		------------------------------------------------
		listcoefficient():

		List the coefficient in Coefficient

		------------------------------------------------
		"""
		m = 0
		label1 = ""
		label2 = ""
		for i in self.coefficients:
			if i != 0:
				print('Z'+str(m+1)+' = ',i,self.zernikelist[m])
				label1 = label1 + 'Z'+str(m+1)+' = '+str(i)+"\n"
				label2 = label2 + 'Z'+str(m+1)+' = '+str(i)+"  "
			m = m + 1
		return [label1,label2]

	def zernikelist(self):
		"""
		------------------------------------------------
		zernikelist():

		List all Zernike Polynomials

		------------------------------------------------
		"""
		m = 1
		for i in self.zernikelist:
			print("Z"+str(m)+":"+i)
			m = m + 1

	def zernikesurface(self, label = True, zlim=[], matrix = False):
		"""
		------------------------------------------------
		zernikesurface(self, label_1 = True):

		Return a 3D Zernike Polynomials surface figure

		label_1: default show label

		------------------------------------------------
		"""
		theta = np.linspace(0, 2*np.pi, 100)
		rho = np.linspace(0, 1, 100)
		[u,r] = np.meshgrid(theta,rho)
		X = r*cos(u)
		Y = r*sin(u)
		Z = interferometer.zernikepolar(self.coefficients,r,u)
		fig = plt.figure(figsize=(12, 8), dpi=80)
		ax = fig.gca(projection='3d')
		surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.RdYlGn,
	        linewidth=0, antialiased=False, alpha = 0.6)

		if zlim == []:
			v = max(abs(Z.max()),abs(Z.min()))
			ax.set_zlim(-v*5, v*5)
			cset = ax.contourf(X, Y, Z, zdir='z', offset=-v*5, cmap=cm.RdYlGn)
		else:
			ax.set_zlim(zlim[0], zlim[1])
			cset = ax.contourf(X, Y, Z, zdir='z', offset=zlim[0], cmap=cm.RdYlGn)

		ax.zaxis.set_major_locator(LinearLocator(10))
		ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
		fig.colorbar(surf, shrink=1, aspect=30)


		p2v = round(tools.peak2valley(Z),5)
		rms1 = round(tools.rms(Z),5)

		label_1 = self.listcoefficient()[0]+"P-V: "+str(p2v)+"\n"+"RMS: "+str(rms1)
		if label == True:
			plt.title('Zernike Polynomials Surface',fontsize=18)
			ax.text2D(0.02, 0.1, label_1, transform=ax.transAxes,fontsize=14)
		else:
			pass
		plt.show()

		if matrix == True:
			return Z
		else:
			pass
	def zernikemap(self, label = True):
		"""
		------------------------------------------------
		zernikemap(self, label_1 = True):

		Return a 2D Zernike Polynomials map figure

		label: default show label

		------------------------------------------------
		"""


		theta = np.linspace(0, 2*np.pi, 400)
		rho = np.linspace(0, 1, 400)
		[u,r] = np.meshgrid(theta,rho)
		X = r*cos(u)
		Y = r*sin(u)
		Z = interferometer.zernikepolar(self.coefficients,r,u)
		fig = plt.figure(figsize=(12, 8), dpi=80)
		ax = fig.gca()
		im = plt.pcolormesh(X, Y, Z, cmap=cm.RdYlGn)

		if label == True:
			plt.title('Zernike Polynomials Surface Heat Map',fontsize=18)
			ax.set_xlabel(self.listcoefficient()[1],fontsize=18)
		plt.colorbar()
		ax.set_aspect('equal', 'datalim')
		plt.show()

	def zernikeline(self):
		"""
		------------------------------------------------
		zernikeline()

		Return a 1D cutoff through x and y axis of a 3D
		Zernike Polynomials surface figure
		------------------------------------------------
		"""

		X = np.linspace(-1, 1, 100)
		Y = np.linspace(-1, 1, 100)
		ZX = interferometer.zernikecartesian(self.coefficients,X,0)
		ZY = interferometer.zernikecartesian(self.coefficients,0,Y)
		fig = plt.figure()
		ax = fig.gca()
		plt.plot(X,ZX)
		plt.plot(Y,ZY)
		plt.grid()
		plt.show()

	def zernikematrix(self,l = 100):
		x = np.linspace(-1, 1, l)
		[X,Y] = np.meshgrid(x,x)
		Z = interferometer.zernikecartesian(self.coefficients,X,Y)
		return Z

	def psfcaculator(self,r=1,lambda_1=632*10**(-9),z=0.1):
		"""
		pupil: Exit pupil diameter
		z: Distance from exit pupil to image plane
		r: pupil radius, in unit of lambda
		"""
		pupil = l1 = 200 # exit pupil sample points
		x = np.linspace(-r, r, l1)
		[X,Y] = np.meshgrid(x,x)
		Z = interferometer.zernikecartesian(self.coefficients,X,Y)
		for i in range(len(Z)):
			for j in range(len(Z)):
				if x[i]**2+x[j]**2>r**2:
					Z[i][j] = 0
		d = 400 # background
		A = np.zeros([d,d])
		A[d//2-l1//2+1:d//2+l1//2+1,d//2-l1//2+1:d//2+l1//2+1] = Z
		axis_1 = d//pupil*r
		fig = plt.figure()
		# ax = fig.gca()
		# plt.imshow(A,extent=[-axis_1,axis_1,-axis_1,axis_1],cmap=cm.RdYlGn)
		# ax.set_xlabel('mm',fontsize=14)
		# plt.colorbar()
		# plt.show()

		abbe = np.exp(-1j*2*np.pi*A)
		for i in range(len(abbe)):
			for j in range(len(abbe)):
				if abbe[i][j]==1:
					abbe[i][j]=0
		PSF = fftshift(fft2(fftshift(abbe)))**2
		PSF = PSF/PSF.max()
		return PSF

	def psf(self,r=1,lambda_1=632*10**(-9),z=0.1):
		"""
		------------------------------------------------
		psf()

		Return the point spread function of a wavefront described by
		Zernike Polynomials
		------------------------------------------------
		Input:

		r: exit pupil radius(mm)

		lambda_1: wavelength(m)

		z: exit pupil to image plane distance(m)

		"""
		print(r,lambda_1,z)
		PSF = self.psfcaculator(r=r,lambda_1=lambda_1,z=z)
		fig = plt.figure(figsize=(9, 6), dpi=80)
		plt.imshow(abs(PSF),cmap=cm.RdYlGn)
		plt.colorbar()
		plt.show()
		return 0

	def otf(self,r=1,lambda_1=632*10**(-9),z=0.1):
		PSF = self.psfcaculator(r=r,lambda_1=lambda_1,z=z)
		OTF = fftshift(fft2(PSF))
		return 0


	def mtf(self,r=1,lambda_1=632*10**(-9),z=0.1,matrix = False):
		"""
		Modulate Transfer function
		"""
		PSF = self.psfcaculator(r=r,lambda_1=lambda_1,z=z)
		MTF = fftshift(fft2(PSF))
		MTF = MTF/MTF.max()
		f0 = r/1000/lambda_1/z/10000   # cutoff frequency?
		fig = plt.figure(figsize=(9, 6), dpi=80)
		plt.imshow(abs(MTF),cmap=cm.bwr)
		plt.colorbar()
		plt.show()
		if matrix == True:
			return MTF
		else:
			return 0

	def ptf(self):
		"""
		Phase transfer function
		"""
		PSF = self.psfcaculator()
		PTF = fftshift(fft2(PSF))
		PTF = np.angle(PTF)
		b = 400
		R = (200)**2
		for i in range(b):
			for j in range(b):
				if (i-b/2)**2+(j-b/2)**2>R:
					PTF[i][j] = 0
		plt.imshow(abs(PTF),cmap=cm.rainbow)
		plt.colorbar()
		plt.show()
		return 0



	def twyman_green(self,lambda_1=632,PR=1):
		interferometer.twyman_green(self,lambda_1=lambda_1,PR=PR)


	def hartmann(self,r=1,R=1):
		M = hartmann.hartmann(self,r=r,R=R)
		return M

	def zernike2seidel(self):
		'''
		Ap is the piston aberration,coefficients Ai represent the
		peak value of the corresponding Seidel aberration term,
		'''
		a = [0]+self.coefficients
		#Piston
		Ap = a[1]-sqrt(3)*a[4]+sqrt(5)*a[11]
		#tilt
		At = 2*sqrt((a[2]-sqrt(8)*a[8])**2+(a[3]-sqrt(8)*a[7])**2)
		Bt = arctan2(a[3]-sqrt(8)*a[7],a[2]-sqrt(8)*a[8])*180/np.pi
		#Astigmatism
		Aa = 2*sqrt(6*(a[5]**2+a[6]**2))
		Ba = 0.5*arctan2(a[5],a[6])*180/np.pi
		#defocus
		Ad = 2*(sqrt(3)*a[4]-3*sqrt(5)*a[11]-Aa)
		#Coma
		Ac = 6*sqrt(2*(a[7]**2+a[8]**2))
		Bc = arctan2(a[7],a[8])*180/np.pi
		#Spherical
		As = 6*sqrt(5)*a[11]
		A = [Ap,At,Bt,Ad,Aa,Ba,Ac,Bc,As]


		seidellist=["Piston",
				 	"Tilt",
				 	"Defocus",
				 	"Astigmatism",
				 	"Coma",
				 	"Spherical"]
		Atable = [[Ap,0.0],[At,Bt],[Ad,0.0],[Aa,Ba],[Ac,Bc],[As,0.0]]
		print("                 Magnitude  Angle (Degrees)")
		print("-------------------------------------------")
		for i in range(len(seidellist)):
			print("| {0:>13s} |  {1:>8s}  | {2:>8s}   |".\
			format(seidellist[i],str(round(Atable[i][0],3)),str(round(Atable[i][1],3))))
		print("-------------------------------------------")
		SeidelCoefficient = seidel2.Coefficient(Atable)
		return SeidelCoefficient
	def removepiston(self):
		"""
		Remove piston, it is just same value for whole aberration map
		"""
		Z = self.coefficients
		Z[0] = 0
		return Z
	def removetilt(self):
		"""
		Remove tilt, it is mainly caused by system tilt, not aberration
		on surface
		"""
		tilt = [2,3]
		Z = self.coefficients
		for i in tilt:
			Z[i-1] = 0
		return Z
	def removecoma(self):
		"""
		Remove coma, most of coma is caused by misalinement
		??? Is high order coma also caused by misalinement ???
		"""
		coma = [7,8,16,17,29,30]
		Z = self.coefficients
		for i in coma:
			Z[i-1] = 0
		return Z




def fitting(Z,n,remain3D=False,remain2D=False,barchart=False,interferogram=False,removepiston=True):
	"""
	------------------------------------------------
	fitting(Z,n)

	Fitting an aberration to several orthonormal Zernike
	polynomials.

	Return: n-th Zernike coefficients for a fitting surface aberration
			Zernike coefficients barchart
			Remaining aberration
			Fiting surface plot
	Input:
	Z: A surface or aberration matrix measure from inteferometer
	   or something else.

	n: How many order of Zernike Polynomials you want to fit

	reamin(default==Flase): show the surface after remove fitting
	aberrations.

	removepiston: if remove piston, default = True
	------------------------------------------------
	"""


	fitlist = []
	l = len(Z)
	x2 = np.linspace(-1, 1, l)
	y2 = np.linspace(-1, 1, l)
	[X2,Y2] = np.meshgrid(x2,y2)
	r = np.sqrt(X2**2 + Y2**2)
	u = np.arctan2(Y2, X2)
	for i in range(n):
		C = [0]*i+[1]+[0]*(37-i-1)
		ZF = interferometer.zernikepolar(C,r,u)
		for i in range(l):
			for j in range(l):
				if x2[i]**2+y2[j]**2>1:
					ZF[i][j]=0
		a = sum(sum(Z*ZF))*2*2/l/l/np.pi
		fitlist.append(round(a,3))


	l1 = len(fitlist)
	fitlist = fitlist+[0]*(37-l1)
	Z_new = Z - interferometer.zernikepolar(fitlist,r,u)
	for i in range(l):
		for j in range(l):
			if x2[i]**2+y2[j]**2>1:
				Z_new[i][j]=0

	#plot bar chart of zernike
	if barchart == True:
		fitlist1 = fitlist[0:n]
		index = np.arange(n)
		fig = plt.figure(figsize=(9, 6), dpi=80)
		xticklist = []
		width = 0.6
		for i in index:
			xticklist.append('Z'+str(i+1))
		barfigure = plt.bar(index, fitlist1, width,color = '#2E9AFE',edgecolor = '#2E9AFE')
		plt.xticks( index+width//2, xticklist )
		plt.xlabel('Zernike Polynomials',fontsize=18)
		plt.ylabel('Coefficient',fontsize=18)
		plt.title('Fitting Zernike Polynomials Coefficient',fontsize=18)

		plt.show()
	else:
		pass


	if remain3D == True:

		fig = plt.figure(figsize=(12, 8), dpi=80)
		ax = fig.gca(projection='3d')
		surf = ax.plot_surface(X2, Y2, Z_new, rstride=1, cstride=1, cmap=cm.RdYlGn,
	        linewidth=0, antialiased=False, alpha = 0.6)
		v = max(abs(Z.max()),abs(Z.min()))
		ax.set_zlim(-v, v)
		ax.zaxis.set_major_locator(LinearLocator(10))
		ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
		cset = ax.contourf(X2, Y2, Z_new, zdir='z', offset=-v, cmap=cm.RdYlGn)
		fig.colorbar(surf, shrink=1, aspect=30)
		plt.title('Remaining Aberration',fontsize=18)
		p2v = round(tools.peak2valley(Z_new),5)
		rms1 = round(tools.rms(Z_new),5)
		label_new = "P-V: "+str(p2v)+"\n"+"RMS: "+str(rms1)
		ax.text2D(0.02, 0.1,label_new, transform=ax.transAxes)
		plt.show()
	else:
		pass

	if remain2D == True:
		fig = plt.figure(figsize=(9, 6), dpi=80)
		ax = fig.gca()
		im = plt.pcolormesh(X2, Y2, Z_new, cmap=cm.RdYlGn)
		plt.colorbar()
		plt.title('Remaining Aberration',fontsize=18)
		ax.set_aspect('equal', 'datalim')
		plt.show()
	else:
		pass

	if interferogram == True:
		zernike_coefficient = Coefficient(fitlist)
		interferometer.twyman_green(zernike_coefficient)
	else:
		pass
	if removepiston == True:
		fitlist[0] = 0
	else:
		pass
	C = Coefficient(fitlist)  #output zernike Coefficient class
	tools.zernikeprint(fitlist)
	return fitlist,C

def transform(zernikes, tx, ty, thetaR, scaling=1.0):
	"""
		------------------------------------------------
		transform(zernikes, scaling, tx, ty, thetaR):

		scale (concentric), translation of center, and rotate your zernike polynomial without re-sampling.
		scaling and translation is performed first and then rotation.

		zernikes: Zernike coefficients in standard ANSI order - a list of floats with arbitrary length. The first index
		 is the pivot, then x-tilt, y-tilt, first astigmatism,...
		tx: Horizontal translation
		ty: Vertical translation
		thetaR: rotation in _degrees_
		scaling: concentring scaling factor.

		return: Zernike coefficients in standard ANSI order - same list length as the input.

		Based on: Linda Lundström and Peter Unsbo, "Transformation of Zernike coefficients: scaled, translated, and
		rotated wavefronts with circular and elliptical pupils," J. Opt. Soc. Am. A 24, 569-577 (2007)
		Translated into Python by: Oskar Truffer
		------------------------------------------------
		"""
	C1 = np.array(zernikes)
	etaS = scaling
	etaT = 2 * np.sqrt(np.square(tx) + np.square(ty)) / 1.0
	thetaT = np.arctan2(ty, tx)
	thetaR = thetaR * np.pi / 180
	jnm = len(C1) - 1
	nmax = int(np.ceil((-3 + np.sqrt(9 + 8 * jnm)) / 2))
	jmax = int(nmax * (nmax + 3) / 2)
	S = np.zeros((jmax + 1, 1))
	S[0:len(C1)] = C1.reshape((len(C1), 1))
	C1 = S
	jmaxTuple = (jmax + 1, jmax + 1)
	P = np.zeros(jmaxTuple)
	N = np.zeros(jmaxTuple)
	R = np.zeros(jmaxTuple)
	CC1 = np.zeros((jmax + 1, 1), dtype=np.complex)
	counter = 0
	for m in xrange(-nmax, nmax + 1):
		for n in xrange(np.abs(m), nmax + 1, 2):
			jnm = (m + n * (n + 2)) / 2
			P[counter, jnm] = 1
			N[counter, counter] = np.sqrt(n + 1)
			for s in xrange(0, int((n - abs(m)) / 2 + 1)):
				R[counter - s, counter] = math.pow((-1), s) * math.factorial(n - s) / (
					math.factorial(s) * math.factorial((n + m) / 2 - s) * math.factorial((n - m) / 2 - s))
			if m < 0:
				CC1[jnm] = complex(C1[(-m + n * (n + 2)) / 2], C1[jnm]) / np.sqrt(2)
			elif m == 0:
				CC1[jnm] = complex(C1[jnm], 0)
			else:
				CC1[jnm] = complex(C1[jnm], -C1[(-m + n * (n + 2)) / 2]) / math.sqrt(2)
			counter += 1

	ETA = None
	for m in xrange(-nmax, nmax + 1):
		for n in xrange(abs(m), nmax + 1, 2):
			transd = _transformMatrix(n, m, jmax, etaS, etaT, thetaT, thetaR)
			if ETA is None:
				ETA = np.matmul(P, transd)
			else:
				ETA = np.hstack([ETA, np.matmul(P, transd)])

	C = np.matmul(np.matmul(np.matmul(np.matmul(np.matmul(np.matmul(inv(P), inv(N)), inv(R)), ETA), R), N), P)
	CC2 = np.matmul(C, CC1)
	C2 = np.zeros((jmax + 1, 1))
	for m in xrange(-nmax, nmax + 1):
		for n in xrange(abs(m), nmax + 1, 2):
			jnm = (m + n * (n + 2)) / 2
			if m < 0:
				C2[jnm] = np.imag(CC2[jnm] - CC2[(-m + n * (n + 2)) / 2]) / np.sqrt(2)
			elif m == 0:
				C2[jnm] = np.real(CC2[jnm])
			else:
				C2[jnm] = np.real(CC2[jnm] + CC2[(-m + n * (n + 2)) / 2]) / np.sqrt(2)

	return C2


def _transformMatrix(n, m, jmax, etaS, etaT, thetaT, thetaR):
	Eta = np.zeros((jmax + 1, 1), dtype=np.complex)
	for p in xrange(0, int((n + m) / 2 + 1)):
		for q in xrange(0, int((n - m) / 2 + 1)):
			nnew = n - p - q
			mnew = m - p + q
			jnm = (mnew + nnew * (nnew + 2)) / 2
			Eta[int(math.floor(jnm))] += _ncr((n + m) / 2, p) * _ncr((n - m) / 2, q) * np.power(etaS,
																										n - p - q) * np.power(
				etaT, p + q) * np.exp(complex(0, ((p - q) * (thetaT - thetaR) + m * thetaR)))
	return Eta


def _ncr(n, r):
	r = min(r, n - r)
	if r == 0: return 1
	numer = reduce(op.mul, xrange(n, n - r, -1))
	denom = reduce(op.mul, xrange(1, r + 1))
	return numer // denom
