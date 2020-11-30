"""PyWENO smooth reconstruction example."""

import numpy as np

import sys
sys.path.insert(0, '/home/mcheikh/python/PyWENO')
import pyweno.weno as weno
import matplotlib.pyplot as plt

f = np.cos
F = np.sin

# Pyweno
x = np.linspace(0.0, 2*np.pi, 21)
a = (F(x[1:]) - F(x[:-1]))/(x[1]-x[0])
l, s, wl = weno.reconstruct(a, 5, 'left', return_smoothness=True,return_weights=True)
r  , wr  = weno.reconstruct(a, 5, 'right',return_weights=True)

# Fortran_Weno
l_m  = np.genfromtxt('output/fortran-smooth-left.txt')
r_m  = np.genfromtxt('output/fortran-smooth-right.txt')
s_m  = np.genfromtxt('output/fortran-smooth-smoothness.txt',)
wl_m  = np.genfromtxt('output/fortran-smooth-weights-left.txt',)
wr_m  = np.genfromtxt('output/fortran-smooth-weights-right.txt',)

print 'smoothness error'
print np.max(abs(s - s_m))
print 'weight left error'
print np.max(abs(wl - wl_m))
print 'weight right error'
print np.max(abs(wr - wr_m))
print 'reconstruction left error'
print np.max(abs(l - l_m))
print 'reconstruction right error'
print np.max(abs(r - r_m))

plt.title('pyweno.weno reconstruction and smoothness indicators')

plt.subplot(2,1,1)

plt.plot(x,   f(x), '-k')
plt.plot(x[:-1], l, 'or')
plt.plot(x[1:],  r, 'ob')
plt.plot(x[:-1], l_m, 'xg')
plt.plot(x[1:],  r_m, 'xm')

plt.ylabel('f')
plt.xlabel('x')
plt.legend(['actual', 'left-pyweno', 'right-pyweno', 'left-fortran', 'right-fortran'], fontsize=8)

plt.subplot(2,1,2)

c = 0.5*(x[1:] + x[:-1])

plt.plot(c, s[:,0], 'or')
plt.plot(c, s[:,1], 'ok')
plt.plot(c, s[:,2], 'ob')

plt.ylabel('smoothness')
plt.xlabel('x')
plt.legend(['r=0', 'r=1', 'r=2'])

plt.savefig('smooth.png', format='png')