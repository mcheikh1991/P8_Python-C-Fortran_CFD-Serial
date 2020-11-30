"""PyWENO smooth reconstruction example."""

import numpy as np
import pyweno.weno   as weno

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

f = np.cos
F = np.sin

x = np.linspace(0.0, 2*np.pi, 21)
a = (F(x[1:]) - F(x[:-1]))/(x[1]-x[0])
l, s, wl = weno.reconstruct(a, 5, 'left', return_smoothness=True,return_weights=True)
r  , wr  = weno.reconstruct(a, 5, 'right',return_weights=True)
l_m  = np.genfromtxt('fortran-smooth-left.txt')
r_m  = np.genfromtxt('fortran-smooth-right.txt')
s_m  = np.genfromtxt('fortran-smooth-smoothness.txt',)
wl_m  = np.genfromtxt('fortran-smooth-weights-left.txt',)
wr_m  = np.genfromtxt('fortran-smooth-weights-right.txt',)

print 'smoothness error'
print np.max(s - s_m)

plt.title('pyweno.weno reconstruction and smoothness indicators')

plt.subplot(2,1,1)

plt.plot(x,   f(x), '-k')
plt.plot(x[:-1], l, 'or')
plt.plot(x[1:],  r, 'ob')
plt.plot(x[:-1], l_m, 'xg')
plt.plot(x[1:],  r_m, 'xm')

plt.ylabel('f')
plt.xlabel('x')
plt.legend(['actual', 'left', 'right'])

plt.subplot(2,1,2)

c = 0.5*(x[1:] + x[:-1])

plt.plot(c, s[:,0], 'or')
plt.plot(c, s[:,1], 'ok')
plt.plot(c, s[:,2], 'ob')

plt.ylabel('smoothness')
plt.xlabel('x')
plt.legend(['r=0', 'r=1', 'r=2'])

plt.savefig('smooth.png', format='png')

import matplotlib.pyplot as plt
import numpy as np
levels = np.linspace(MN,MX,LV) # level of colors
CS = plt.contourf(X_Phi, Y_Phi, Phi,levels)
cbar = plt.colorbar(CS)
plt.xlabel(Xlabel)
plt.ylabel(Ylable)
plt.xlim(Xlim)
plt.ylim(Ylim)