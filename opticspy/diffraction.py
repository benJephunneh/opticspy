import numpy as np
from . import tools as tools

def fresnel(aperture,z = 2,lambda1 = 660*10**(-9)):
	aperturelist = ['rectangle','doublerectangle','frame','doublecircle','ring']
	if aperture.type == 'circle':
		n = aperture.background
		d = aperture.d
		D = aperture.d*aperture.scale
		scale = aperture.scale
		Nf = D**2/(4*lambda1*z);  #Fresnel number. Smaller is better for single-DFT Fresnel
		print("Fresnel number = ", Nf)
	elif aperture.type in aperturelist:
		print(aperture.type)
		n = aperture.background
		scale = aperture.scale
	else:
		print("No this kind aperture for fresnel diffraction")

	x1 = np.linspace(-n/2+1,n/2,n)
	[x,y] = np.meshgrid(x1,x1)
	# Single-DFT
	e1 = np.exp(1j*2*np.pi/lambda1*(x**2+y**2)/2/z*((scale)**2))
	diffraction = np.fft.fftshift(np.fft.fft2(aperture.aper*e1))

	extent = n*scale
	tools.apershow(diffraction, extent = extent)
	return diffraction

def fraunhofer(aperture, z = 2, lambda1 = 660*10**(-9)):
	"""
	Fraunhofer diffraction
	"""
	diffraction = 1j*np.exp(1j*2*np.pi/lambda1*z)/lambda1/z*np.fft.fftshift(np.fft.fft2(aperture.aper))

	extent = aperture.background*aperture.scale
	tools.apershow(diffraction, extent)
	return diffraction