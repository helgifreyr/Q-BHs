from data import *
from pylab import *
import sys
import grapefruit

rc('font', family='serif', size=8)

def plot_gogn(r,g,b):
    gogn = get_data('1.0','*','*','*')
    B_gogn_1st = get_Bstar_data('1.0','*','1st')
    B_gogn_2nd = get_Bstar_data('1.0','*','2nd')
    e_gogn_1st = get_extremal_data('1.0','*','1st','*')
    e_gogn_2nd = get_extremal_data('1.0','*','2nd','*')
    ws = gogn[:,18]
    Ms = gogn[:,7]
    e_ws_1st = e_gogn_1st[:,18]
    e_Ms_1st = e_gogn_1st[:,7]
    e_ws_2nd = e_gogn_2nd[:,18]
    e_Ms_2nd = e_gogn_2nd[:,7]
    B_ws_1st = B_gogn_1st[:,1]
    B_Ms_1st = B_gogn_1st[:,5]
    B_ws_2nd = B_gogn_2nd[:,1]
    B_Ms_2nd = B_gogn_2nd[:,5]
    B_ws = concatenate((B_ws_1st[::-1],B_ws_2nd))
    B_Ms = concatenate((B_Ms_1st[::-1],B_Ms_2nd))
    e_ws = concatenate((e_ws_1st,e_ws_2nd[::-1]))
    e_Ms = concatenate((e_Ms_1st,e_Ms_2nd[::-1]))
    colB = grapefruit.Color.NewFromRgb(r, g, b)
    col = colB.MonochromeScheme()[0]
    colE = colB.MonochromeScheme()[1]
    plot(ws,Ms,'.',c=col.html,ms=1.5,label=r'Q-BHs')
    plot(e_ws,e_Ms,'-',c=colE.html,ms=1.5,label=r'eQ-BHs')
    plot(B_ws,B_Ms,'-',c=colB.html,ms=2.5,label=r'Q-balls')

plot_gogn(1,0,0)


xlabel(r'$w$')
ylabel(r'$M$')
xlim(0.75,1.00)
# curve = genfromtxt('zero-mode')
# curve_w = curve[:,0]
# curve_M = curve[:,1]
# plot(curve_w,curve_M,'m')
# no_int = genfromtxt('no-interaction-extremal')
# no_int_w = no_int[:,0]
# no_int_M = no_int[:,1]
# plot(no_int_w,no_int_M,'k')
W = linspace(0.75,1.0,100)
M = lambda w: 1/(2*w)
plot(W,M(W),'k')
legend()
savefig('w-vs-M-all.png')
