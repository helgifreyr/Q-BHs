from data import get_data
from pylab import *
import sys
rc('font', family='serif', size=8)

m = sys.argv[1].replace('x','*')
w = sys.argv[2].replace('x','*')
branch = sys.argv[3].replace('x','*')
rh = sys.argv[4].replace('x','*')
c2 = sys.argv[5].replace('x','*')

gogn = get_data(m,w,branch,rh,c2)
ws = gogn[:,18]
rhs = gogn[:,0]
Js = gogn[:,5]
THs = gogn[:,6]
Mcs = gogn[:,7]
AHs = gogn[:,8]
Mints = gogn[:,15]
suptitle(r'$m='+m+',w='+w+',rh='+rh+',c2='+c2+'$')
style='.'
subplot(321)
xlabel(r'$r_H$')
ylabel(r'$J$')
plot(rhs,-Js,style,ms=1.0)
subplot(322)
xlabel(r'$r_H$')
ylabel(r'$T_H$')
ylim(min(THs),sum(THs)/len(THs)/2)
plot(rhs,THs,style,ms=1.0)
subplot(323)
xlabel(r'$r_H$')
ylabel(r'$M$')
plot(rhs,Mcs,style,ms=1.0)
subplot(324)
xlabel(r'$r_H$')
ylabel(r'$A_H$')
plot(rhs,AHs,style,ms=1.0)
subplot(325)
xlabel(r'$r_H$')
ylabel(r'$\frac{M_{int}}{M}$')
plot(rhs,-Mints/Mcs,style,ms=1.5)
# xlim(0,0.31)
savefig('m='+m+'-w='+w+'-rh='+rh+'-c2='+c2+'.png')
