import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.insert(0, '/home/mcheikh/python/PyWENO')
import pyweno.weno as weno
import matplotlib.pyplot as plt

def calculate_face_value(phi,phi_bl,phi_br,phi_bl_type,phi_br_type):

	k = WENO_INFO['k'] # k-th order for WENO
	type_weno = WENO_INFO['type'] # 'JS' or 'M' or 'Z'

	nx, ny = phi
	n      = [nx,ny]
	
	phi_r  = np.zeros((nx+1, ny))
	phi_l  = np.zeros((nx+1, ny))

	# Finding the Face values
	j = 0
	for i in range(ny):
		# Loop over the x direction (y is fixed):
		phi_1D = phi[:,i]    

		phi_boundary1 = phi_bl[i]
		phi_boundary2 = phi_br[i]

		phi_boundary1_type = phi_bl_type[i]
		phi_boundary2_type = phi_br_type[i]



		#if j==1 and  i == 51:
		#	print phi_1D

		# Apply Boundary Conditions:
		if phi_boundary1_type == 'Periodic':
			phi_1D =  np.concatenate((phi_1D[-(k+1)/2:],phi_1D) , axis = 0)
		elif phi_boundary1_type == 'Dirichlet':
			phi_1D =  np.concatenate((phi_boundary1*(np.ones((k+1)/2)),phi_1D) , axis=0)
		elif phi_boundary1_type == 'ZeroGradient':
			phi_1D =  np.concatenate( (np.ones((k+1)/2)*phi_1D[0],phi_1D) , axis=0)
		else:
			print 'Undefined Boundary Condition'
			raise

		#if j==1 and  i == 51:
		#	print phi_1D

		if phi_boundary2_type == 'Periodic':
			phi_1D =  np.concatenate((phi_1D,phi_1D[:(k+1)/2]), axis = 0)
		elif phi_boundary2_type == 'Dirichlet':
			phi_1D =  np.concatenate((phi_1D,phi_boundary2*(np.ones((k+1)/2))) , axis=0)
		elif phi_boundary1_type == 'ZeroGradient':
			phi_1D =  np.concatenate((phi_1D,np.ones((k+1)/2)*phi_1D[-1]),axis=0)
		else:
			print 'Undefined Boundary Condition'
			raise

		# Finding the face value
		phi_l_1D = weno.reconstruct(phi_1D, k, 'left' ,weno_type=type_weno)[(k+1)/2:-(k+1)/2]
		phi_r_1D = weno.reconstruct(phi_1D, k, 'right',weno_type=type_weno)[(k+1)/2:-(k+1)/2]


		#if j==1 and  i == 51:
		#	print weno.reconstruct(phi_1D, k, 'left' ,weno_type=type_weno)
		#	print '-----'

		# Fixing the face value at the boundaries
		if phi_boundary1_type == 'Periodic':
			phi_r_1D = np.append(phi_r_1D[-1],phi_r_1D)
		elif phi_boundary1_type == 'Dirichlet':
			phi_r_1D = np.append(phi_boundary1,phi_r_1D)
		elif phi_boundary1_type == 'ZeroGradient':
			phi_r_1D = np.append(phi_1D[0],phi_r_1D)
			phi_l_1D[0] = phi_1D[0]
		else:
			print 'Undefined Boundary Condition'
			raise

		if phi_boundary2_type == 'Periodic':
			phi_l_1D = np.append(phi_l_1D ,phi_l_1D[0])
		elif phi_boundary2_type == 'Dirichlet':
			phi_l_1D = np.append(phi_l_1D ,phi_boundary2)
		elif phi_boundary2_type == 'ZeroGradient':
			phi_l_1D = np.append(phi_l_1D ,phi_1D[-1])
			phi_r_1D[-1] = phi_1D[-1]
		else:
			print 'Undefined Boundary Condition'
			raise

		# Saving the Face Value
		phi_l[i,:] = phi_l_1D
		phi_r[i,:] = phi_r_1D

	return phi_l , phi_r

nx,ny = 100,50
x = np.linspace(0,np.pi*2,nx)
y = np.linspace(0,np.pi,ny)

phi = np.zeros([nx,ny])
xm = np.zeros([nx,ny])
ym = np.zeros([nx,ny])
for i in range(len(x)):
	for j in range(len(y)):
		phi[i,j] = np.cos(x[i])*np.sin(y[j])
		xm[i,j]  = x[i]
		ym[i,j]  = y[j]

phi_bl = np.zeros(ny)
phi_br = np.zeros(ny)
phi_bl_type = ['Periodic' for i in range(ny)]
phi_br_type = ['Periodic' for i in range(ny)]
phi_l , phi_r =  calculate_face_value(phi,phi_bl,phi_br,phi_bl_type,phi_br_type)

plt.contourf(xm,ym,phi)
plt.show()