import numpy as np
from math import *
#from pycse import bvp
import matplotlib.pyplot as plt
from numpy import *

timestop = 0.00125
u = 1.0
ut = u*timestop
npoints = zeros(6)
error = zeros(6)

for i in arange(0,6,1):
	npoints[i] = 16*2**i

def ComputeAccuracy(n):
	line = loadtxt('LinePlot%03d.txt'%n)
	x   = line[:,0]
	rho = line[:,1]
	n = size(line[:,0])
	errorval = 0.0
	for i in arange(0,n,1):
        	rhoexact = 0.01*exp(-100.0*((x[i]-ut)-0.5)**2)
        	errorval = errorval + (rho[i]-rhoexact)**2
	errorval = (errorval/n)**0.5
	print errorval
	return errorval

for i in arange(0,6,1):
	n = npoints[i]
	error[i] = ComputeAccuracy(n)

dx = zeros(6)

dx[0] = 1.0/16
dx[1] = 1.0/32
dx[2] = 1.0/64
dx[3] = 1.0/128
dx[4] = 1.0/256
dx[5] = 1.0/512

plt.loglog(1/dx[:],error[:],'or',markersize=10)

slope = zeros([2,2])

slope[0,0] = 10**1.3
slope[0,1] = 10**(-5.0)
slope[1,0] = 10**2.3
slope[1,1] = 10**(-7.3)

plt.plot(slope[:,0],slope[:,1],'k')

plt.title('Order of accuracy screwup',fontsize=30)


plt.show()

