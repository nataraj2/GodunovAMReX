import numpy as np
from math import *
#from pycse import bvp
import matplotlib.pyplot as plt
from numpy import *

timestop = 0.125
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
        	#errorval = errorval + abs(rho[i]-rhoexact)
	errorval = (errorval/n)**0.5
	#errorval = sum(errorval)/n
	print errorval
	return errorval

for i in arange(0,6,1):
	n = npoints[i]
	error[i] = ComputeAccuracy(n)

dx = zeros(6)

for i in arange(0,6,1):
	dx[i] = 1.0/(16*2**i)

plt.loglog(1/dx[:],error[:],'or',markersize=10)

slope = zeros([2,2])

order = log(error[1]/error[2])/log(dx[1]/dx[2])
print "Order = ", order

startval = log10(error[1])
print error[2], startval, 10**(startval)
endval = startval - order
print endval

slope[0,0] = 10**(log10(1/dx[1]))
slope[0,1] = 10**(startval)
slope[1,0] = 10**(1+log10(1/dx[1]))
slope[1,1] = 10**(endval)


plt.plot(slope[:,0],slope[:,1],'k')

plt.title('Order of accuracy = %f'%order,fontsize=30)
#figname = './Images/Order_FullScheme.png'
figname = './Images/Order_NoLimAndMono.png'
#figname = './Images/Order_OnlyLimNoMono.png'
#figname = './Images/Order_NoLimOnlyMono.png'
plt.savefig(figname)

plt.show()

