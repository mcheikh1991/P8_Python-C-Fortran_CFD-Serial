import numpy as np
import sys
import os
import pyweno.weno as weno

class Field_Variable(object):
	"""
	A class that takes a 2D array to solve a 2D problem
	"""

	def __init__(self,name='Variable'):
		self.name = name
		self.num_interior_elements = 0
		self.num_boundary_elements = 0
		pass

	def __str__(self) :
		""" Return string representation of Ray.
		"""
		return self.name

	def set_interior_value(self, interior_array):
		if isinstance(interior_array, np.ndarray):
			self.interior_value = interior_array
			self.num_interior_elements = interior_array.size
		else:
			print("Interior should be of type numpy.ndarray")
			raise
		pass


	def set_boundary_value(self, boundary_array, boundary_type, location):
		"""
		setup the boundary value and type for each boundary cell, if periodic 
		then the value will not matter
		"""

		All_Locations = ['left','right','up','down','bottom','top']
		location = location.lower() 

		if All_Locations.count(location) == 0:
			print("Unknown Location, should be one of the following", All_Locations)
			raise
		else:
			# Incase of a uniform boundary
			if  type(boundary_type) == str:
				boundary_type = [boundary_type for i in range(len(boundary_array))]
			
			if location == 'left':
				self.boundary_left_value = boundary_array
				self.boundary_left_type  = boundary_type
			elif location == 'right':
				self.boundary_right_value = boundary_array
				self.boundary_right_type  = boundary_type
			elif location == 'up' or location == 'top':
				self.boundary_up_value = boundary_array
				self.boundary_up_type  = boundary_type
			elif location == 'down' or location == 'bottom':
				self.boundary_down_value = boundary_array
				self.boundary_down_type  = boundary_type
		pass

	def __add__(self,other):
		if isinstance(self, other.__class__):
			A = self.interior_value + other.interior_value
		else:
			A = self.interior_value + other
		return A  

	def __sub__(self,other):
		if isinstance(self, other.__class__):
			A = self.interior_value - other.interior_value
		else:
			A = self.interior_value - other
		return A  

	def __mul__(self,other):
		if isinstance(self, other.__class__):
			A = self.interior_value * other.interior_value
		else:
			A = self.interior_value * other
		return A

	def __truediv__(self,other):
		if isinstance(self, other.__class__):
			A = self.interior_value / other.interior_value
		else:
			A = self.interior_value / other
		return A

	def __pow__(self,other):
		if isinstance(self, other.__class__):
			A = self.interior_value ** other.interior_value
		else:
			A = self.interior_value ** other
		return A

	def copy(self):
		from copy import deepcopy
		return deepcopy(self)

	def print_to_file(self,filename):
		f = open(filename,'w')
		for j in range(len(self.interior_value)):
			for i in range(len(self.interior_value[j,:])):
				f.write(str(self.interior_value[j,i]))
				f.write(',')
			f.write('\n')
		f.close()

	def set_face_value(self,other,location=''):
		if location   == 'right':
			self.face_right_value = other
		elif location == 'left':
			self.face_left_value  = other
		elif location == 'up'or location == 'top':
			self.face_up_value    = other
		elif location == 'down'or location == 'bottom':
			self.face_down_value  = other
		else:
			print "Unknown location given"
			raise
		pass


	def return_rl_faces(self):
		return np.array([self.face_right_value,self.face_left_value])


	def return_ud_faces(self):
		return np.array([self.face_up_value,self.face_down_value])


	def calculate_face_value(self,WENO_INFO,direction=''):

		if ['x'].count(direction.lower()) == 1:
			# Loop over the x direction (y is fixed):
			Looping_Direction = [0]
		elif ['y'].count(direction.lower()) == 1:
			# Loop over the y direction (x is fixed):
			Looping_Direction = [1]
		else :
			# Loop over the both directions:
			Looping_Direction = [0,1]


		k = WENO_INFO['k'] # k-th order for WENO
		type_weno = WENO_INFO['type'] # 'JS' or 'M' or 'Z'

		nx, ny = self.interior_value.shape
		n      = [nx,ny]
		
		phi_r  = np.zeros((ny, nx+1))
		phi_l  = np.zeros((ny, nx+1))
		phi_u  = np.zeros((ny+1, nx))
		phi_d  = np.zeros((ny+1, nx))


		# Finding the Face values
		for j in Looping_Direction:
			for i in range(n[j]):
				if   j == 0 : 
					# Loop over the x direction (y is fixed):
					phi_1D =  self.interior_value[i,:]    

					phi_boundary1 = self.boundary_left_value[i]
					phi_boundary2 = self.boundary_right_value[i]

					phi_boundary1_type = self.boundary_left_type[i]
					phi_boundary2_type = self.boundary_right_type[i]

				elif j == 1 : 
					# Loop over the y direction (x is fixed):
					phi_1D =  self.interior_value[:,i]    

					phi_boundary1 = self.boundary_down_value[i]
					phi_boundary2 = self.boundary_up_value[i]

					phi_boundary1_type = self.boundary_down_type[i]
					phi_boundary2_type = self.boundary_up_type[i]


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
				if   j == 0 : # Loop over the x direction (y is fixed):
					phi_l[i,:] = phi_l_1D
					phi_r[i,:] = phi_r_1D
				elif j == 1 : # Loop over the y direction (x is fixed):
					phi_d[:,i] = phi_l_1D
					phi_u[:,i] = phi_r_1D

				#print self.name 
				#if self.name == 'y-velocity' and i== 51:
				#	print 	phi_d[:,i] 
				#	print   phi_u[:,i]
				#	print '----------'


			if   j == 0 : # Loop over the x direction (y is fixed):
				self.face_right_value = phi_r
				self.face_left_value  = phi_l
			elif j == 1 : # Loop over the y direction (x is fixed):
				self.face_up_value    = phi_u
				self.face_down_value  = phi_d

		pass


