import numpy as np
from numpy import sin as sin
from numpy import cos as cos
import matplotlib.pyplot as plt

def doubleslit(b=0.1,a=0.4,lambda_1=632,z=0.5):
    """
    Return a Young's doubleslit(Frauhofer Diffraction)
    Input:
    --------------------------------
    b: slit of width in mm
    a: slit separation of in mm
    lambda_1: wavelength in nm
    z: slit-to-screen distance in m.
    """
    lambda_1 = float(lambda_1)
    theta = np.linspace(-0.04,0.04,1000)
    theta1 = np.ones(100)
    [theta,theta1] = np.meshgrid(theta,theta1)
    beta = np.pi*(b/1000)/(lambda_1/(10**9))*sin(theta)
    alpha = np.pi*(a/1000)/(lambda_1/(10**9))*sin(theta)
    y = 4*(sin(beta)**2/(beta**2)*cos(alpha)**2)
    fig = plt.figure(1,figsize=(12,8), dpi=80)
    plt.imshow(-y)
    plt.set_cmap('Greys')
    plt.show()
    
    theta = np.linspace(-0.04,0.04,1000)
    beta = np.pi*(b/1000)/(lambda_1/(10**9))*sin(theta)
    alpha = np.pi*(a/1000)/(lambda_1/(10**9))*sin(theta)
    y = 4*(sin(beta)**2/(beta**2)*cos(alpha)**2)
    y1 = 4*sin(beta)**2/(beta**2)
    fig = plt.figure(2,figsize=(12, 8), dpi=80)
    plt.plot(theta*z*1000,y)
    plt.plot(theta*z*1000,y1,"g--")
    plt.show()
