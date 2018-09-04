import numpy as np
import matplotlib.pyplot as plt
def third(s1,s2,s3,s4,s5):

	"""
	Third order aberrations:
	Ray aberrations
	Field curve
	Distortion

	input: third order aberration coefficient
		   sigma 1~5

	output: third order aberration graph
	"""

	print("third order aberration")
	py = np.linspace(-1,1,100)
	px = np.linspace(0,1,50)

	height = [1,0.7,0]
	count = 0
	ax = []
	maxTan = 0
	maxSag = 0


	fig = plt.figure(1)
	for h in height:
		Tan = s1*py**3+3*s2*h*py**2+(3*s3+s4)*h**2.*py+s5*h**3
		ax.append(plt.subplot2grid((3, 3), (count, 0), colspan=2))
		plt.plot(py, Tan)
		if maxTan < max(abs(Tan)): maxTan = max(abs(Tan))
		if count == 0: plt.title('TANGENTIAL')
		plt.axis([-1, 1, -maxTan, maxTan])
		if count == len(height)-1: plt.xlabel('\n' + r'$\rho_y$',fontsize=20)
		plt.ylabel('h = '+str(h),fontsize=15)
		plt.grid(True)

		Sag = s1*px**3+(s3+s4)*h**2*px
		ax.append(plt.subplot2grid((3, 3), (count, 2)))
		plt.plot(px, Sag)
		if maxSag < max(abs(Sag)): maxSag = max(abs(Sag))
		plt.axis([0, 1, -maxSag, maxSag])
		if count == 0: plt.title('SAGITTAL')
		if count == len(height)-1: plt.xlabel('\n' + r'$\rho_x$',fontsize=20)
		plt.grid(True)

		count = count + 1

	fig.set_tight_layout(True)
	plt.show()

def fieldcurve(sigma3 = 0.05, sigma4 = -0.05, FNO = 10, H = 20):
	"""
	sigma3  Astigmatism Coefficient
	sigma4  Petzval Coefficient
	FNO     F-number
	H       Image Height
	"""
	uak = -1.00/(2*FNO)   # maginal ray angle
	h = np.linspace(0,1,40)
	XP = -sigma4/uak*h**2
	XT = -(3*sigma3+sigma4)/uak*h**2
	XS = -(sigma3+sigma4)/uak*h**2
	fig = plt.figure(figsize=(6, 8), dpi=80)
	plt.plot(XP, h*H, 'b-*', label='P')
	plt.plot(XT, h*H, 'b--', label='T')
	plt.plot(XS, h*H, 'b', label='S')
	plt.xlabel('Surface sag(mm)',fontsize=18)
	plt.ylabel('Real image height(mm)',fontsize=18)
	legend = plt.legend(loc='lower left', shadow=True, fontsize='x-large')
	plt.title(r'$\sigma3 = $'+str(round(sigma3,4))+' '+r'$\sigma4 = $'+str(sigma4),fontsize=18)
	#plt.axis([-16, 5, 0, H])
	plt.grid(b=True, which='both', color='0.65',linestyle='--')
	plt.show()
	return 0