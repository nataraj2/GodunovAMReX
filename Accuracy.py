import numpy as np
from math import *
#from pycse import bvp
import matplotlib.pyplot as plt
from numpy import *

timestop = loadtxt('timestop.txt')
u = 1.0
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
	ut = u*timestop[i]
	error[i] = ComputeAccuracy(n)

dx = zeros(6)

for i in arange(0,6,1):
	dx[i] = 1.0/(16*2**i)

plt.loglog(1/dx[:],error[:],'or',markersize=10)

slope = zeros([2,2])

startat = 0
endat = startat+1

order =1.8#log(error[startat]/error[endat])/log(dx[startat]/dx[endat])
print "Order = ", order

startval = log10(error[startat])
print error[2], startval, 10**(startval)
endval = startval - order
print endval

slope[0,0] = 10**(log10(1/dx[startat]))
slope[0,1] = 10**(startval)
slope[1,0] = 10**(1+log10(1/dx[startat]))
slope[1,1] = 10**(endval)


plt.plot(slope[:,0],slope[:,1],'k')

cfl = loadtxt('cfl.txt')

plt.title('CFL=%f, Order = %.02f'%(cfl,order),fontsize=25)
plt.xlabel('No. of points',fontsize=20)
plt.ylabel('L2 norm',fontsize=20)
figname = './Images/OrderForCFL_%f.png'%cfl
#figname = './Images/Order_NoLimAndNoMono.png'
#figname = './Images/Order_OnlyLimNoMono.png'
#figname = './Images/Order_NoLimOnlyMono.png'
plt.savefig(figname)

plt.show()

