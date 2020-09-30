import numpy as np
from math import *
#from pycse import bvp
import matplotlib.pyplot as plt
from numpy import *

xmin = 0.0
xmax = 1.0
cfl = 0.25
u = 1.0
dt = 5e-4;#cfl*dx/u
niter = 250

def vanLeer(a,b,c):
    small_qty_sq = (1.e-10)*(1.e-10);

    dsc = 0.5e0*(b - c);
    dsl = 2.e0*(a - c);
    dsr = 2.e0*(b - a);
    zero = 0.0
      
    if(dsl*dsr > small_qty_sq):
	if(dsc>0):
		return min(abs(dsc),min(abs(dsl),abs(dsr)))
	else:
		return -1.0*min(abs(dsc),min(abs(dsl),abs(dsr)))
    else:
		return zero;

def getfpandfm(s,i,dx):
	sm2 = s[i-2]
	sm1 = s[i-1]
	s0 =  s[i]
	sp1 = s[i+1]
	sp2 = s[i+2]

        d1 = vanLeer(s0,sp1,sm1);
        d2 = vanLeer(sm1,s0,sm2);

	sedge1 = 0.5*(s0+sm1) - 1.0/6.0*(d1-d2)
	sedge1 = min(max(sedge1,min(s0, sm1)),max(s0,sm1));

	d1 = vanLeer(sp1,sp2,s0);
    	d2 = vanLeer(s0,sp1,sm1);

        sedge2 = 0.5*(sp1 + s0) - 1.0/6.0*(d1-d2)
        sedge2 = min(max(sedge2,min(s0, sp1)),max(s0,sp1));

	sm = sedge1
	sp = sedge2

	#if ((sedge2-s0)*(s0-sedge1) < 0.e0):
        #	sp = s0;
        #	sm = s0;
    	#elif (abs(sedge2-s[i]) >= 2.0*abs(sedge1-s0)):
      	#	sp = 3.0*s0 - 2.0*sedge1;
	#elif (abs(sedge1-s[i]) >= 2.0*abs(sedge2-s0)):
	#	sm = 3.0*s0 - 2.0*sedge2;
	
	s6 = 6.0*s0 - 3.0*(sm + sp);

        sigmap = u*dt/dx;
        sigmam = u*dt/dx;

        flux = sp - (0.5*sigmap)*((sp - sm) - (1.e0 -2.e0/3.e0*sigmap)*s6)

	return flux

def AdvectionEquation(n):
	x = zeros(n)
	s = zeros(n)
	sold = zeros(n)
	snew = zeros(n)
	fL = zeros(n+1)
	dx = (xmax-xmin)/n
	for i in arange(0,n,1):
		dx = (xmax-xmin)/n
                x[i] = xmin + (i+0.5)*dx
                s[i] = 0.01*exp(-100.0*(x[i]-0.5)**2)
                snew[i] = s[i]
                sold[i] = s[i]
		fL[i] = 0.0
	for iteration in arange(0,niter+1,1):
		print iteration, iteration*dt
		for i in arange(2,n-3,1):
			fL[i+1] = getfpandfm(sold,i,dx)	
	
		for i in arange(2,n-3,1):
			s[i] = sold[i] + u*dt/dx*(fL[i]-fL[i+1])
		
		for i in arange(2,n-3,1):
			sold[i] = s[i]	
	
		if(mod(iteration,10)==0):
			plt.plot(x,s)
			#figname = './Images/AdvectionGodunov%04d.png'%(iteration/10)
			#plt.savefig(figname)
			plt.pause(0.001)
			plt.clf()

	filename = 'LinePlot%03d.txt'%(n)		
	f = open(filename,"w")
	for i in arange(0,n,1):
		f.write(("%.16f %.16f\n")%(x[i],s[i]))
	f.close()

for i in arange(0,6,1):
	n = 16*2**i
	AdvectionEquation(n)

