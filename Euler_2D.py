"""PyWENO smooth reconstruction example."""

import numpy as np
from classes.Class_Field_Variable import *
from functions.Func_Initialize import *
from functions.Func_Thermo import *
from functions.Func_Plotting import *
import fortran_functions.Euler_Flux_Schemes_2D as EF # Fortran Code compiled using: f2py -c twod_euler_fluxes_v2.f90 -m Euler_Flux_Schemes_2D
EF_2D = EF.convective_schemes_2d

import fortran_functions.Euler_Flux_Schemes_1D as EF # Fortran Code compiled using: f2py -c oned_euler_fluxes_v5.f90 -m Euler_Flux_Schemes_1D
EF_1D = EF.convective_schemes_1d
'''
We are going to solve the 
Euler equations here:
	d(rho)/dt + d(rho*u)/dx + d(rho*u)/dy 
	du/dt + d(rho u^2 + p)/dx + d(rho*u*v)/dy = 0 
	dv/dt + d(rho*u*v)/dx + d(rho v^2 + p)/dy = 0 
	dE/dt + d(u(E+p))/dx + d(v(E+p))/dy = 0 
	p = (gamma - 1) * (E - rho/2 * (u^2 + v^2))
'''
def print_file(Case,Name):
    f = open(Name,'w')
    for i in range(len(Case)):
        f.write(str(Case[i]))
        f.write('\n')
    f.close()

autoplot = 10
dt = 0.0002*10 #1e-5
nx = 101
ny = nx
nt = int(0.5 / dt)
scheme_cons = 2
scheme_time = 1

# WENO Informatiom
k = 5
Weno_Info = {}
Weno_Info['k']     = k # Order of WENO
Weno_Info['type']  = 'JS' # JS or M or Z

# Transport Properties
gamma = 1.4 # Specific heat ratio
CP = 2.5
CV = 1.786
Field_Constants ={}
Field_Constants['gamma']= gamma
Field_Constants['CP']= gamma
Field_Constants['CV']= gamma

#Domain and Time
x    = np.linspace(-1.0,1.0,nx)
y    = np.linspace(-1.0,1.0,ny)

dx = x[1] - x[0]
dy = y[1] - y[0]
X,Y = np.meshgrid(x, y)


plottinglimits = [-0.1,1.1,0.8,3.5]

Mesh_Info = {}
Mesh_Info['X'] = X
Mesh_Info['Y'] = Y
Mesh_Info['y'] = y
Mesh_Info['x'] = x
Mesh_Info['ny'] = ny
Mesh_Info['nx'] = nx

u   = Field_Variable('x-velocity')
v   = Field_Variable('y-velocity')
rho = Field_Variable('density')
p   = Field_Variable('pressure')
E   = Field_Variable('Total_Energy')
e   = Field_Variable('Internal_Energy')

# Inital and Boundary Conditions
u,v,rho,p,E,e = Init_Bound_Variables2(u,v,rho,p,E,e,Mesh_Info,Field_Constants)


# KNP Approach

rhoU  = Field_Variable('x-density-velocity')
rhoV  = Field_Variable('y-density-velocity')

rhoU.interior_value = rho*u
rhoV.interior_value = rho*v

rhoU.set_boundary_value( rho.boundary_right_value*u.boundary_right_value, 'ZeroGradient', 'right')
rhoU.set_boundary_value( rho.boundary_left_value*u.boundary_left_value, 'ZeroGradient', 'left')
rhoU.set_boundary_value( rho.boundary_up_value*u.boundary_up_value, 'ZeroGradient', 'up')
rhoU.set_boundary_value( rho.boundary_down_value*u.boundary_down_value, 'ZeroGradient', 'down')

rhoV.set_boundary_value( rho.boundary_right_value*v.boundary_right_value, 'ZeroGradient', 'right')
rhoV.set_boundary_value( rho.boundary_left_value*v.boundary_left_value, 'ZeroGradient', 'left')
rhoV.set_boundary_value( rho.boundary_up_value*v.boundary_up_value, 'ZeroGradient', 'up')
rhoV.set_boundary_value( rho.boundary_down_value*v.boundary_down_value, 'ZeroGradient', 'down')


rho.calculate_face_value(Weno_Info)

# x-direction Flux
x_phi   = np.zeros(rho.face_right_value.shape)
x_phiUp = np.zeros(rho.face_right_value.shape)
x_phiVp = np.zeros(rho.face_right_value.shape)
x_phiEp = np.zeros(rho.face_right_value.shape)

# y-direction Flux
y_phi   = np.zeros(rho.face_up_value.shape)
y_phiUp = np.zeros(rho.face_up_value.shape)
y_phiVp = np.zeros(rho.face_up_value.shape)
y_phiEp = np.zeros(rho.face_up_value.shape)

plot_contour(X,Y,rho.interior_value,0,'Results/Cont-Density-%i'%0,[0.1,1.1])
plot_contour3D(X,Y,rho.interior_value,0,'Results/Cont-3D-Density-%i'%0,[0.1,1.1])


rho_1 	= rho.copy()
u_1		= u.copy()
v_1		= v.copy()
E_1		= E.copy()
p_1		= p.copy()
e_1		= e.copy()
rhoU_1  = rhoU.copy()
rhoV_1  = rhoV.copy()

print_file(rho.interior_value[50,:],'rho_%d' % 0)
print_file(u.interior_value[50,:],'u_%d' % 0)
print_file(v.interior_value[50,:],'v_%d' % 0)
print_file(E.interior_value[50,:],'E_%d' % 0)
print_file(p.interior_value[50,:],'p_%d' % 0)

for t in range(1,nt+1):
	print 'time: ', t*dt

	# First Part
	#--------------
	rho.calculate_face_value(Weno_Info)
	u.calculate_face_value(Weno_Info)
	v.calculate_face_value(Weno_Info)
	p.calculate_face_value(Weno_Info)

	rhoU.set_interior_value(rho*u)
	rhoV.set_interior_value(rho*v)
	E.set_interior_value( rho*0.5*(u**2.0 + v**2.0) + p*((gamma-1.0)**(-1)) ) 

	rhoU.set_face_value(rho.face_right_value*u.face_right_value,'right')
	rhoU.set_face_value(rho.face_left_value*u.face_left_value,'left')
	rhoU.set_face_value(rho.face_up_value*u.face_up_value,'up')
	rhoU.set_face_value(rho.face_down_value*u.face_down_value,'down')

	rhoV.set_face_value(rho.face_right_value*v.face_right_value,'right')
	rhoV.set_face_value(rho.face_left_value*v.face_left_value,'left')
	rhoV.set_face_value(rho.face_up_value*v.face_up_value,'up')
	rhoV.set_face_value(rho.face_down_value*v.face_down_value,'down')


	E.set_face_value(rho.face_right_value*0.5*(u.face_right_value**2.0 + \
		v.face_right_value**2.0) + p.face_right_value/(gamma - 1.0) ,'right')
	E.set_face_value(rho.face_left_value*0.5*(u.face_left_value**2.0 + \
		v.face_left_value**2.0) + p.face_left_value/(gamma - 1.0),'left')
	E.set_face_value(rho.face_up_value*0.5*(u.face_up_value**2.0 + \
		v.face_up_value**2.0) + p.face_up_value/(gamma - 1.0),'up')
	E.set_face_value(rho.face_down_value*0.5*(u.face_down_value**2.0 + \
		v.face_down_value**2.0) + p.face_down_value/(gamma - 1.0),'down')

	
	# Loop Over the faces in x-direction (y-fixed)
	for j in range(len(rho.face_right_value)):
		for i in range(len(rho.face_right_value[j])):
			RIGHT = [rho.face_right_value[j,i],rhoU.face_right_value[j,i],rhoV.face_right_value[j,i],E.face_right_value[j,i]]
			LEFT  = [rho.face_left_value[j,i],rhoU.face_left_value[j,i],rhoV.face_left_value[j,i],E.face_left_value[j,i]]
			[x_phi[j,i],x_phiUp[j,i],x_phiVp[j,i],x_phiEp[j,i]] = EF_2D.roe(RIGHT,LEFT,1,0)

	# Loop Over the faces in y-direction (x-fixed)
	for j in range(len(rho.face_up_value)):
		for i in range(len(rho.face_up_value[j])):
			RIGHT = [rho.face_up_value[j,i],rhoU.face_up_value[j,i],rhoV.face_up_value[j,i],E.face_up_value[j,i]]
			LEFT  = [rho.face_down_value[j,i],rhoU.face_down_value[j,i],rhoV.face_down_value[j,i],E.face_down_value[j,i]]
			[y_phi[j,i],y_phiUp[j,i],y_phiVp[j,i],y_phiEp[j,i]] = EF_2D.roe(RIGHT,LEFT,0,1)

	L_rho   = (-1/(dx))*(x_phi[:,1:]   - x_phi[:,:-1])   + (-1/(dy))*(y_phi[1:,:] - y_phi[:-1,:])
	L_rho_u = (-1/(dx))*(x_phiUp[:,1:] - x_phiUp[:,:-1]) + (-1/(dy))*(y_phiUp[1:,:] - y_phiUp[:-1,:])
	L_rho_v = (-1/(dx))*(x_phiVp[:,1:] - x_phiVp[:,:-1]) + (-1/(dy))*(y_phiVp[1:,:] - y_phiVp[:-1,:])
	L_E     = (-1/(dx))*(x_phiEp[:,1:] - x_phiEp[:,:-1]) + (-1/(dy))*(y_phiEp[1:,:] - y_phiEp[:-1,:])

	rho_1.set_interior_value(rho.interior_value + dt*L_rho)
	u_1.set_interior_value((rho.interior_value*u.interior_value + dt*L_rho_u)/rho_1.interior_value)
	v_1.set_interior_value((rho.interior_value*v.interior_value + dt*L_rho_v)/rho_1.interior_value)
	E_1.set_interior_value(E.interior_value + dt*L_E)
	p_1.set_interior_value(p_idealgas(E_1.interior_value,rho_1.interior_value,u_1.interior_value,v_1.interior_value,gamma))
	e_1.set_interior_value(e_internal_energy(E_1.interior_value,rho_1.interior_value,u_1.interior_value,v_1.interior_value))

	# Second Part
	#--------------
	rhoU_1.set_interior_value(rho_1*u_1)
	rhoV_1.set_interior_value(rho_1*v_1)

	rho_1.calculate_face_value(Weno_Info)
	u_1.calculate_face_value(Weno_Info)
	v_1.calculate_face_value(Weno_Info)	
	p_1.calculate_face_value(Weno_Info)

	rhoU_1.set_face_value(rho_1.face_right_value*u_1.face_right_value,'right')
	rhoU_1.set_face_value(rho_1.face_left_value*u_1.face_left_value,'left')
	rhoU_1.set_face_value(rho_1.face_up_value*u_1.face_up_value,'up')
	rhoU_1.set_face_value(rho_1.face_down_value*u_1.face_down_value,'down')

	rhoV_1.set_face_value(rho_1.face_right_value*v_1.face_right_value,'right')
	rhoV_1.set_face_value(rho_1.face_left_value*v_1.face_left_value,'left')
	rhoV_1.set_face_value(rho_1.face_up_value*v_1.face_up_value,'up')
	rhoV_1.set_face_value(rho_1.face_down_value*v_1.face_down_value,'down')

	
	E_1.set_face_value(rho_1.face_right_value*0.5*(u_1.face_right_value**2.0 + \
		v_1.face_right_value**2.0) + p_1.face_right_value/(gamma - 1.0) ,'right')
	E_1.set_face_value(rho_1.face_left_value*0.5*(u_1.face_left_value**2.0 + \
		v_1.face_left_value**2.0) + p_1.face_left_value/(gamma - 1.0),'left')
	E_1.set_face_value(rho_1.face_up_value*0.5*(u_1.face_up_value**2.0 + \
		v_1.face_up_value**2.0) + p_1.face_up_value/(gamma - 1.0),'up')
	E_1.set_face_value(rho_1.face_down_value*0.5*(u_1.face_down_value**2.0 + \
		v_1.face_down_value**2.0) + p_1.face_down_value/(gamma - 1.0),'down')

	# Loop Over the faces in x-direction (y-fixed)
	for j in range(len(rho_1.face_right_value)):
		for i in range(len(rho_1.face_right_value[j])):
			RIGHT = [rho_1.face_right_value[j,i],rhoU_1.face_right_value[j,i],rhoV_1.face_right_value[j,i],E_1.face_right_value[j,i]]
			LEFT  = [rho_1.face_left_value[j,i],rhoU_1.face_left_value[j,i],rhoV_1.face_left_value[j,i],E_1.face_left_value[j,i]]
			[x_phi[j,i],x_phiUp[j,i],x_phiVp[j,i],x_phiEp[j,i]] = EF_2D.roe(RIGHT,LEFT,1,0)
			
	# Loop Over the faces in y-direction (x-fixed)
	for j in range(len(rho_1.face_up_value)):
		for i in range(len(rho_1.face_up_value[j])):
			RIGHT = [rho_1.face_up_value[j,i],rhoU_1.face_up_value[j,i],rhoV_1.face_up_value[j,i],E_1.face_up_value[j,i]]
			LEFT  = [rho_1.face_down_value[j,i],rhoU_1.face_down_value[j,i],rhoV_1.face_down_value[j,i],E_1.face_down_value[j,i]]
			[y_phi[j,i],y_phiUp[j,i],y_phiVp[j,i],y_phiEp[j,i]] = EF_2D.roe(RIGHT,LEFT,0,1)

	
	L_rho_1   = (-1/(dx))*(x_phi[:,1:]   - x_phi[:,:-1])   + (-1/(dy))*(y_phi[1:,:] - y_phi[:-1,:])
	L_rho_u_1 = (-1/(dx))*(x_phiUp[:,1:] - x_phiUp[:,:-1]) + (-1/(dy))*(y_phiUp[1:,:] - y_phiUp[:-1,:])
	L_rho_v_1 = (-1/(dx))*(x_phiVp[:,1:] - x_phiVp[:,:-1]) + (-1/(dy))*(y_phiVp[1:,:] - y_phiVp[:-1,:])
	L_E_1     = (-1/(dx))*(x_phiEp[:,1:] - x_phiEp[:,:-1]) + (-1/(dy))*(y_phiEp[1:,:] - y_phiEp[:-1,:])


	rho_new = 0.5*rho.interior_value + 0.5*rho_1.interior_value+ 0.5*dt*L_rho_1
	u_new   = (0.5*u.interior_value*rho.interior_value + 0.5*u_1.interior_value*rho_1.interior_value +  0.5*dt*L_rho_u_1)/rho_new
	v_new   = (0.5*v.interior_value*rho.interior_value + 0.5*v_1.interior_value*rho_1.interior_value +  0.5*dt*L_rho_v_1)/rho_new
	E_new   = 0.5*E.interior_value + 0.5*E_1.interior_value + 0.5*dt*L_E_1 
	p_new   = p_idealgas(E_new,rho_new,u_new,v_new,gamma)
	e_new   = e_internal_energy(E_new,rho_new,u_new,v_new)

	rho.set_interior_value(rho_new)
	u.set_interior_value(u_new)
	v.set_interior_value(v_new)
	E.set_interior_value(E_new)
	p.set_interior_value(p_new)
	e.set_interior_value(e_new)

	if t % autoplot == 0:
		plot_contour(X,Y,rho.interior_value,t,'Results/Cont-Density-%i'%t,[0.1,1.1])
		plot_contour3D(X,Y,rho.interior_value,t,'Results/Cont-3D-Density-%i'%t,[0.1,1.1])
		plot_contour3D(X,Y,u.interior_value,t,'Results/Cont-3D-Vel-U-%i'%t,[-0.1,1.1])
		plot_contour3D(X,Y,v.interior_value,t,'Results/Cont-3D-Vel-V-%i'%t,[-0.1,1.1])
		plot_fig(x,[rho.interior_value[2,:],rho.interior_value[51,:]],int(t/autoplot),k,['Near-The-Top','Center'],[-1.1,1.1,0.05,1.1],'Results/Density',ylabel_name='$\\rho$',Leg_loc='upper right',marker=False)
		plot_fig(x,[u.interior_value[2,:],u.interior_value[51,:]],int(t/autoplot),k,['Near-The-Top','Center'],[-1.1,1.1,-0.1,1.1],'Results/X-velocity',ylabel_name='$u_x$',Leg_loc='upper right',marker=False)
		plot_fig(x,[v.interior_value[2,:],v.interior_value[51,:]],int(t/autoplot),k,['Near-The-Top','Center'],[-1.1,1.1,-0.1,1.1],'Results/Y-velocity',ylabel_name='$u_y$',Leg_loc='upper right',marker=False)
		plot_fig(x,[p.interior_value[2,:],p.interior_value[51,:]],int(t/autoplot),k,['Near-The-Top','Center'],[-1.1,1.1,0.05,1.1],'Results/Pressure',ylabel_name='$p$',Leg_loc='upper right',marker=False)
		
	if t*dt == 0.2:
		Exact_Sol = np.genfromtxt('Results/Exact_Solution_1D_SOD/output-0.2-N_1000',skip_header=2)
		cell_ex, x_ex,  rho_ex, p_ex, u_ex = Exact_Sol[:,0] , Exact_Sol[:,1] , Exact_Sol[:,2], Exact_Sol[:,3], Exact_Sol[:,4]
		
		plot_fig(x,[rho.interior_value[2,:],rho.interior_value[51,:]],int(t/autoplot),k,['Near-The-Top','Center'],[-1.1,1.1,0.05,1.1],'Results/Density',ylabel_name='$\\rho$',x_opt=x_ex,y_opt=rho_ex,name_opt='Analytical Sol',Leg_loc='upper right',marker=False)
		plot_fig(x,[rho.interior_value[2,:],rho.interior_value[51,:]],int(t/autoplot),k,['Near-The-Top','Center'],[-0.3,-0.20,0.95,1.02],'Results/DensityA',ylabel_name='$\\rho$',x_opt=x_ex,y_opt=rho_ex,name_opt='Analytical Sol',Leg_loc='lower right')
		plot_fig(x,[rho.interior_value[2,:],rho.interior_value[51,:]],int(t/autoplot),k,['Near-The-Top','Center'],[-0.1,0.2,0.35,0.52],'Results/DensityB',ylabel_name='$\\rho$',x_opt=x_ex,y_opt=rho_ex,name_opt='Analytical Sol',Leg_loc='lower right')
		plot_fig(x,[rho.interior_value[2,:],rho.interior_value[51,:]],int(t/autoplot),k,['Near-The-Top','Center'],[0.1,0.28,0.25,0.46],'Results/DensityC',ylabel_name='$\\rho$',x_opt=x_ex,y_opt=rho_ex,name_opt='Analytical Sol',Leg_loc='lower right')
		