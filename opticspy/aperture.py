import numpy as np

from . import diffraction as diffraction
from . import tools as tools

class Aperture():
	def __init__(self, background):
		self.background = background

	def show(self):

		"""
		Show aperture figure

		Output
		------------------------------------
		Aperture figure

		"""

		print("---------show aperture--------")
		extent = self.scale * self.background
		tools.apershow(self.aper, extent)

	def fresnel(self,z = 2,lambda1 = 660*10**(-9)):
		"""
		Compute the fresnel diffraction pattern.

		Output
		------------------------------------
		Diffraction pattern figure

		"""
		print("---------Fresnel Diffraction----------")
		diffraction.fresnel(self,z,lambda1)
		return 0
	def fraunhofer(self,z = 2,lambda1 = 660*10**(-9)):
		"""
		Compute the frauhofer diffraction pattern.

		Output
		------------------------------------
		Diffraction pattern figure

		"""
		print("---------Fraunhofer Diffraction----------")
		diffraction.fraunhofer(self,z,lambda1)
		return 0

	def otf(self):
		"""
		Compute an aperture's otf
		"""
		print("-------------OTF---------------")
		aperfft = np.fft.fftshift(np.fft.fft2(self.aper))**2
		aper_OTF = np.fft.fftshift(np.fft.fft2(aperfft))
		tools.apershow(aper_OTF,extent = 0)
		return 0


class Circle(Aperture):
	"""
	Build a circle aperture example

	Example
	------------------------------------
	aperture = opticspy.aper.Circle(200,50)

 	Parameters
	------------------------------------
	background: int
				Square background

	d: int
				aperture pixel diameter
	D: int
				aperture real diameter

	"""
	def __init__(self, background=500, d=200, D = 0.01, scale = 0.01/200):
		self.type = 'circle'
		self.background = n = background
		self.d = d
		self.D = d*scale
		self.scale = scale
		radius = d/2
		self.aper = np.zeros([n,n])
		aper1 = tools.circle_aperture(d)
		self.aper[(n//2-d//2):(n//2-d//2+d),(n//2-d//2):(n//2-d//2+d)] = aper1

class DoubleCircle(Aperture):
	def __init__(self, background=500, d=50, D=0.01, separation = 100, scale=0.01/200):
		self.type = 'doublecircle'
		self.background = n = background
		self.d = DoubleRectangle
		self.D = D
		self.scale = scale
		self.separation = s = separation
		radius = d/2
		self.aper = np.zeros([n,n])

		aper1 = tools.circle_aperture(d)
		self.aper[(n//2-d//2):(n//2-d//2+d),(n//2-s//2-d//2):(n//2-s//2+d//2)] = aper1
		self.aper[(n//2-d//2):(n//2-d//2+d),(n//2+s//2-d//2):(n//2+s//2+d//2)] = aper1

class Ring(Aperture):
	def __init__(self, background=500, outside=200, inside=100, scale=0.01/200):
		self.type = 'ring'
		self.background = n = background
		self.outside = outside
		self.inside = inside
		self.scale = scale
		self.aper = np.zeros([n,n])

		aper1 = tools.circle_aperture(inside)
		aper2 = tools.circle_aperture(outside)

		aper2[(outside//2-inside//2):(outside//2-inside//2+inside),(outside//2-inside//2):(outside//2-inside//2+inside)] = -1*(aper1-1)
		self.aper[(n//2-outside//2):(n//2-outside//2+outside),(n//2-outside//2):(n//2-outside//2+outside)] = aper2

class Rectangle(Aperture):
	def __init__(self, background=500, height=200, width=200, scale=0.01/200):
		"""
		Build a rectangle aperture instance

		Example
		-----------
		aperture = opticspy.aper.Rectangle(200,20,40)

	 	Parameters
		-----------
		background: int
					Square background
		height: int
					aperture height
		width: int
					aperture width
		"""
		self.type = 'rectangle'
		n = self.background = background
		self.height = height
		self.width = width
		self.scale = scale
		#matrix_1 = [height,width]
		aper1 = np.ones([height,width])
		self.aper = np.zeros([n,n])
		self.aper[(n//2-height//2):(n//2-height//2+height),(n//2-width//2):(n//2-width//2+width)] = aper1

class DoubleRectangle(Aperture):
	def __init__(self, background=500, height=50, width=2, separation=4, scale=0.01/200):
		"""
		Build a DoubleRectangle aperture instance, could use as a doubleslit aperture
		"""
		self.type = "doublerectangle"
		n = self.background = background
		self.height = height
		self.width = width
		self.separation = separation
		self.scale = scale
		self.aper = np.zeros([n,n])
		aper1 = np.ones([height,width])
		self.aper[(n//2-height//2):(n//2-height//2+height),(n//2-width//2-separation//2):(n//2-width//2+width-separation//2)] = aper1
		self.aper[(n//2-height//2):(n//2-height//2+height),(n//2-width//2+separation//2):(n//2-width//2+width+separation//2)] = aper1

class Frame(Aperture):
	def __init__(self, background=500, outside=200, inside=100, scale=0.01/200):
		self.type = "frame"
		n = self.background = background
		self.outside = outside
		self.inside = inside
		self.scale = scale
		self.aper = np.zeros([n,n])
		aper1 = np.ones([outside,outside])
		aper2 = np.zeros([inside,inside])
		self.aper[(n//2-outside//2):(n//2-outside//2+outside),(n//2-outside//2):(n//2-outside//2+outside)] = aper1
		self.aper[(n//2-inside//2):(n//2-inside//2+inside),(n//2-inside//2):(n//2-inside//2+inside)] = aper2
