"""PyWENO smooth reconstruction example."""

import numpy as np
import matplotlib.pyplot as plt

problem = 1 # 1 for smooth and 2 for discontinuous
use_pyweno = 1

for problem in [1,2]:

	if problem == 1:
		f = np.cos
		F = np.sin

		x = np.linspace(0.0, 2*np.pi, 21)
		a = (F(x[1:]) - F(x[:-1]))/(x[1]-x[0]) # a = d (sin(x)) / dx 	

	elif problem == 2:
		def f(x):
			r = np.zeros(x.shape)
			i = x > 0
			r[i] = np.cos(x[i])
			i = x <= 0
			r[i] = np.sin(x[i])
			return r

		def F(x):
			r = np.zeros(x.shape)
			i =x > 0
			r[i] = np.sin(x[i])
			i = x <= 0
			r[i] = -np.cos(x[i])
			return r
		x = np.linspace(-2*np.pi, 2*np.pi, 41)
		a = (F(x[1:]) - F(x[:-1]))/(x[1]-x[0])
		a[20:21] = (F(x[21:22]) - 0.0)/(x[1]-x[0]) # fix middle cell average

	if use_pyweno == 1:
		import pyweno.weno
		l, sl       = pyweno.weno.reconstruct(a, 5, 'left', return_smoothness=True)
		r, sr , wr  = pyweno.weno.reconstruct(a, 5, 'right', return_smoothness=True , return_weights=True)
		m, sm , wm  = pyweno.weno.reconstruct(a, 5, 'middle', return_smoothness=True , return_weights=True)
			
	'''
	Everything below is an expansion of what is occuring in pyweno.weno.reconstruct(a, 5, 'right' and 'left')
	---------------------------------------------------------------
	'''

	q = a			 		   # cell-averaged unknown to reconstruct
	k = 5			 		   # order of reconstruction (odd values from 5 to 11)
	points = 'right'   		   # reconstruction points at the left edge of each cell
	return_smoothness = False  # use smoothness indicators *smoothness* (optional)
	smoothness = None		   # Default
	weights	= None		
	return_weights = False,
	squeeze=True
	N = q.shape[0]	# 20
	k = (k+1)/2	   # 3

	# Updating points if it was gauss
	if points == 'gauss':
		points = 'gauss_legendre'

	# validate points and n
	if points in [ 'left', 'right', 'middle' ]:
		n = 1
	elif n is None:
		n = k

	if points in [ '+/-' ]:
		n = 2		   # n = 1 because the points is right
		
	# Part 1 (smoothness):
	# --------------------
		
	smoothness = np.zeros((N,k)) # np.zeros([20,3])
	# since k = 3 we will go to weno\_smoothness003.c which is defined below

	#sigma = np.zeros([3])
	for i in range(3,N-3):
		
		sigma0 =  \
			(10/3.0) * q[i + 0] * q[i + 0] + \
			(-31/3.0)* q[i + 0] * q[i + 1] + \
			(11/3.0) * q[i + 0] * q[i + 2] + \
			(25/3.0) * q[i + 1] * q[i + 1] + \
			(-19/3.0)* q[i + 1] * q[i + 2] + \
			(4/3.0)  * q[i + 2] * q[i + 2]
			
		sigma1 = \
			(4/3.0)   * q[i - 1] * q[i - 1] + \
			(-13/3.0) * q[i - 1] * q[i + 0] + \
			(10.0/6.0)* q[i - 1] * q[i + 1] + \
		  	(13/3.0)  * q[i + 0] * q[i + 0] + \
			(-13/3.0) * q[i + 0] * q[i + 1] + \
			(4/3.0)   * q[i + 1] * q[i + 1]
				
		sigma2 = \
			(4/3.0)  * q[i - 2] * q[i - 2] + \
			(-19/3.0)* q[i - 2] * q[i - 1] + \
			(11/3.0) * q[i - 2] * q[i + 0] + \
			(25/3.0) * q[i - 1] * q[i - 1] + \
			(-31/3.0)* q[i - 1] * q[i + 0] + \
			(10/3.0) * q[i + 0] * q[i + 0]
		
		smoothness[i] = [sigma0,sigma1,sigma2]
			
	# Part 2 (weights):
	# --------------------

	weights_right = np.zeros((N,k))
	weights_left  = np.zeros((N,k))
	weights_middle = np.zeros((N,k,2))

	for i in range(3,N-3):

		sigma0 = smoothness[i,0]
		sigma1 = smoothness[i,1]
		sigma2 = smoothness[i,2]
		
		# Calculating the weights for right (v_{i+1/2})
		acc = 0.0;
		omega0 = (+0.3) / ((sigma0 + 1.0e-6) * (sigma0 + 1.0e-6));
		acc = acc + omega0;
		omega1 = (+0.6) / ((sigma1 + 1.0e-6) * (sigma1 + 1.0e-6));
		acc = acc + omega1;
		omega2 = (+0.1) / ((sigma2 + 1.0e-6) * (sigma2 + 1.0e-6));
		acc = acc + omega2;
		# Normalize 
		omega0 = (omega0) / (acc); 
		omega1 = (omega1) / (acc);
		omega2 = (omega2) / (acc);
		weights_right[i , 0 ] = omega0;
		weights_right[i , 1 ] = omega1;
		weights_right[i , 2 ] = omega2;
		
		# Calculating the weights for left (v_{i-1/2})
		acc = 0.0;
		omega0 = (+0.1) / ((sigma0 + 1.0e-6) * (sigma0 + 1.0e-6));
		acc = acc + omega0;
		omega1 = (+0.6) / ((sigma1 + 1.0e-6) * (sigma1 + 1.0e-6));
		acc = acc + omega1;
		omega2 = (+0.3) / ((sigma2 + 1.0e-6) * (sigma2 + 1.0e-6));
		acc = acc + omega2;
		# Normalize 
		omega0 = (omega0) / (acc); 
		omega1 = (omega1) / (acc);
		omega2 = (omega2) / (acc);
		weights_left[i , 0 ] = omega0;
		weights_left[i , 1 ] = omega1;
		weights_left[i , 2 ] = omega2;

		# Calculating the weights for middle (v_{i}) [Beed to find source for it]
		acc = 0.0;
		omega0p = ((+0.1125) / (+2.675)) / ((sigma0 + 1.0e-6) * (sigma0 + 1.0e-6));
		acc = acc + omega0p;
		omega1p = ((+2.45) / (+2.675)) / ((sigma1 + 1.0e-6) * (sigma1 + 1.0e-6));
		acc = acc + omega1p;
		omega2p = ((+0.1125) / (+2.675)) / ((sigma2 + 1.0e-6) * (sigma2 + 1.0e-6));
		acc = acc + omega2p;
		# Normalize 
		omega0p = (omega0p) / (acc);
		omega1p = (omega1p) / (acc);
		omega2p = (omega2p) / (acc);

		acc = 0.0;
		omega0m = ((+0.225) / (+1.675)) / ((sigma0 + 1.0e-6) * (sigma0 + 1.0e-6));
		acc = acc + omega0m;
		omega1m = ((+1.225) / (+1.675)) / ((sigma1 + 1.0e-6) * (sigma1 + 1.0e-6));
		acc = acc + omega1m;
		omega2m = ((+0.225) / (+1.675)) / ((sigma2 + 1.0e-6) * (sigma2 + 1.0e-6));		
		# Normalize 
		acc = acc + omega2m;
		omega0m = (omega0m) / (acc);
		omega1m = (omega1m) / (acc);
		omega2m = (omega2m) / (acc);

		weights_middle[i , 0 , 0] = omega0p;
		weights_middle[i , 0 , 1] = omega0m;
		weights_middle[i , 1 , 0] = omega1p;
		weights_middle[i , 1 , 1] = omega1m;
		weights_middle[i , 2 , 0] = omega2p;
		weights_middle[i , 2 , 1] = omega2m;
		  
	# Part 3 (reconstruct):
	# --------------------

	q_right = np.zeros((N,n))
	q_left  = np.zeros((N,n))
	q_middle  = np.zeros((N,n))

	for i in range(3,N-3):
		# Reconstruction for the right (v_{i+1/2})
		omega0 = weights_right[i,0]
		omega1 = weights_right[i,1]
		omega2 = weights_right[i,2]
		
		qr0 = (1/3.0) * q[i + 0] + (5/6.0)* q[i + 1] + (-1/6.0) * q[i + 2]
		qr1 = (-1/6.0) * q[i - 1] + (5/6.0)* q[i + 0] + (1/3.0) * q[i + 1]
		qr2 = (1/3.0) * q[i - 2] + (-7/6.0)* q[i - 1] + (11/6.0) * q[i + 0]	  			
		
		q_right[i] = (omega0) * (qr0) + (omega1) * (qr1) + (omega2) * (qr2)
		
		# Reconstruction for the left (v_{i-1/2})
		omega0 = weights_left[i,0]
		omega1 = weights_left[i,1]
		omega2 = weights_left[i,2]
		 
		qr0 = (11/6.0) * q[i + 0] + (-7/6.0)* q[i + 1] + (1/3.0) * q[i + 2]
		qr1 = (1/3.0) * q[i - 1] + (5/6.0)* q[i + 0] + (-1/6.0)  * q[i + 1]
		qr2 = (-1/6.0) * q[i - 2] + (5/6.0)* q[i - 1] + (1/3.0)  * q[i + 0]
		
		q_left[i] = (omega0) * (qr0) + (omega1) * (qr1) + (omega2) * (qr2)

		# Reconstruction for the middle (v_{i})
		omega0p = weights_middle[i , 0 , 0] 
		omega0m = weights_middle[i , 0 , 1] 
		omega1p = weights_middle[i , 1 , 0]
		omega1m = weights_middle[i , 1 , 1]
		omega2p = weights_middle[i , 2 , 0]
		omega2m = weights_middle[i , 2 , 1]

		fr0 = (23/24.0)*q[i + 0] + (1/12.0)*q[i + 1 ] + (-1/24.0) * q[i + 2]
		fr1 = (-1/24.0)*q[i - 1] + (13/12.0)*q[i + 0] + (-1/24.0) * q[i + 1]
		fr2 = (-1/24.0)*q[i - 2] + (1/12.0)*q[i - 1]  + (23/24.0) * q[i + 0]
		q_middle[i] = ((+2.675) * ((omega0p) * (fr0) + (omega1p) * (fr1) + (omega2p) * (fr2))) - \
			((+1.675) * ((omega0m) * (fr0) + (omega1m) * (fr1) + (omega2m) * (fr2)));
	# Plot Data Below Below
		
	x2 = np.linspace(x[0], x[-1], 1001)
	plt.figure(problem)
	plt.plot(x2,f(x2),'C0',label='correct solution')
	#plt.plot(x,f(x),'8b')
	plt.plot(x[:-1], q_left, 'sC1' ,label='WENO left k=5')  # Similiar to pyweno
	plt.plot(x[1:],  q_right, 'oC2',label='WENO right k=5') # Similiar to pyweno
	plt.plot((x[1:]+x[:-1])/2.0,  q_middle, '*C5',label='WENO middle k=5') # Similiar to pyweno
	if use_pyweno == 1:
		plt.plot(x[:-1], l, 'xC3' ,label='pyWENO left k=5')  # Similiar to pyweno
		plt.plot(x[1:],  r, 'dC4',label='pyWENO right k=5') # Similiar to pyweno
		plt.plot( (x[1:]+x[:-1])/2.0,  m, '+C6',label='pyWENO middle k=5') # Similiar to pyweno
	plt.ylabel('f(x)')
	plt.legend()
	plt.savefig('test_%d.png'%problem)


