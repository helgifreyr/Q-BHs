from data import *
from pylab import *
import sys
rc('font', family='serif', size=8)

m = sys.argv[1].replace('x','*')
w = sys.argv[2].replace('x','*')
branch = sys.argv[3].replace('x','*')
rh = sys.argv[4].replace('x','*')
c2 = sys.argv[5].replace('x','*')

gogn = get_data(m,w,branch,rh,c2)
B_gogn = get_Bstar_data(m,w,branch,c2)
ws = gogn[:,18]
Mcs = gogn[:,7]
B_ws = B_gogn[:,1]
B_Mcs = B_gogn[:,5]
title(r'$m='+m+',w='+w+',rh='+rh+',c2='+c2+'$')
xlabel(r'$w$')
ylabel(r'$M$')
xlim(0.75,1.0)
plot(ws,Mcs,'b.',ms=1.5)
plot(B_ws,B_Mcs,'r.',ms=2.5)
curve = genfromtxt('zero-mode')
curve_w = curve[:,0]
curve_M = curve[:,1]
plot(curve_w,curve_M,'m')
W = linspace(0.75,1.0,100)
M = lambda w: 1/(2*w)
plot(W,M(W),'k')
savefig('w-vs-M-'+c2+'.png')
