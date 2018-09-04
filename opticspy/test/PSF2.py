import numpy as np
from numpy import sqrt as sqrt
from numpy import cos as cos
from numpy import sin as sin
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from matplotlib.ticker import LinearLocator as LinearLocator
from matplotlib.ticker import FormatStrFormatter as FormatStrFormatter
from numpy.fft import fftshift as fftshift
from numpy.fft import ifftshift as ifftshift
from numpy.fft import fft2 as fft2

l1 = 100
#Generate test surface matrix from a detector
x = np.linspace(-1, 1, l1)
y = np.linspace(-1, 1, l1)
[X,Y] = np.meshgrid(x,y)
r = sqrt(X**2+Y**2)
#Z = 20*X
Z = sqrt(14)*(8*X**4-8*X**2*r**2+r**4)*(6*r**2-5)
for i in range(len(Z)):
	for j in range(len(Z)):
		if x[i]**2+y[j]**2>1:
			Z[i][j]=0

d = 400
A = np.zeros([d,d])
A[d/2-49:d/2+51,d/2-49:d/2+51] = Z
plt.imshow(A)
plt.show()
def exp_func(a):
	if a == 0:
		return 0
	else:
		return np.exp(1j*2*np.pi*a)
exp_func1 = np.vectorize(exp_func,otypes=[np.complex64])
abbe = exp_func1(A)

fig = plt.figure(2)
AP = abs(fftshift(fft2(fftshift(abbe))))**2
AP = AP/AP.max()
plt.imshow(AP)
plt.show()