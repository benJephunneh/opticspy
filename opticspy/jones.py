"""Simple Jones Vector class
see http://en.wikipedia.org/wiki/Jones_calculus
currently does not normalize the output
"""
from __future__ import division as division
import numpy as np
from numpy import pi as pi


def rotator(angle):
    s = np.sin(angle)
    c = np.cos(angle)
    return np.matrix([[c,-s],[s,c]])

    
class Component(np.matrix):

    def new(self,matrix):
        return super(Component,self).new(self,matrix)
    
    def rotate(self, angle):
        return rotator(angle)*self*rotator(-angle)

    
def HalfWavePlate():
    return Component([[1,0],[0,-1]])
def QuaterWavePlate():
    return Component([[1,0],[0,1j]])
def Birefringence( w1, w2):
    return Component([[np.exp(1j*w1),0],[0,np.exp(1j*w2)]])


def Hpol():
    return np.matrix([[1],[0]])
def Vpol():
    return np.matrix([[0],[1]])
def D1pol():
    return np.matrix([[1],[1]])/np.sqrt(2)
def D2pol():
    return np.matrix([[1],[-1]])/np.sqrt(2)
def C1pol():
    return np.matrix([[1],[1j]])/np.sqrt(2)
def C2pol():
    return np.matrix([[1],[-1j]])/np.sqrt(2)
    

def PolarizerH():
    return Component(Hpol()*Hpol().T)
def PolarizerV():
    return Component(Vpol()*Vpol().T)


if name == "__main__":
    #Usage example
    print(QuaterWavePlate().rotate(pi/4)*HalfWavePlate().rotate(pi/8)*D1pol())

