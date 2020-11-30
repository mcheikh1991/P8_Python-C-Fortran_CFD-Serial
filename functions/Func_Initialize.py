import numpy as np

def Init_Bound_Variables(u,v,rho,p,E,e,Mesh_Info,Field_Constants):
	'''
	This function initialises the variables and sets the boundary condition
	'''
	x  = Mesh_Info['x']
	X  = Mesh_Info['X']

	y  = Mesh_Info['y']
	Y  = Mesh_Info['Y']

	nx = Mesh_Info['nx']
	ny = Mesh_Info['ny']

	gamma = Field_Constants['gamma']

	# Inital Conditions
	p0   = np.array([ 1.0   ,     0.4,  0.0439,    0.15])
	rho0 = np.array([ 0.2579,  0.5197,  0.1072,       1])
	u0   = np.array([ 0.0   , -0.7259, -0.7259,     0.0])
	v0   = np.array([ 0.0   ,     0.0, -1.4045, -1.4045])
	E0   = p0/(gamma-1) + (rho0/2.0) * (u0**2.0 + v0**2.0)	
	e0   = E0/rho0 -  (u0**2.0 + v0**2.0)

	p1   = p0[0]*(((X>=0.5)*1 + (Y>=0.5)*1)==2) + p0[1]*(((X<0.5)*1  + (Y>=0.5)*1)==2) +  \
	     + p0[2]*(((X<0.5)*1  + (Y<0.5)*1)==2)  + p0[3]*(((X>=0.5)*1 + (Y<0.5)*1)==2)

	rho1 = rho0[0]*(((X>=0.5)*1 + (Y>=0.5)*1)==2) + rho0[1]*(((X<0.5)*1  + (Y>=0.5)*1)==2) +  \
	     + rho0[2]*(((X<0.5)*1  + (Y<0.5)*1)==2)  + rho0[3]*(((X>=0.5)*1 + (Y<0.5)*1)==2)

	u1   = u0[0]*(((X>=0.5)*1 + (Y>=0.5)*1)==2) + u0[1]*(((X<0.5)*1  + (Y>=0.5)*1)==2) +  \
	     + u0[2]*(((X<0.5)*1  + (Y<0.5)*1)==2)  + u0[3]*(((X>=0.5)*1 + (Y<0.5)*1)==2)

	v1   = v0[0]*(((X>=0.5)*1 + (Y>=0.5)*1)==2) + v0[1]*(((X<0.5)*1  + (Y>=0.5)*1)==2) +  \
	     + v0[2]*(((X<0.5)*1  + (Y<0.5)*1)==2)  + v0[3]*(((X>=0.5)*1 + (Y<0.5)*1)==2)

	E1   = p1/(gamma-1) + (rho1/2.0) * (u1**2.0 + v1**2.0)
	e1   = E1/rho1 -  (u1**2.0 + v1**2.0)

	u.set_interior_value(u1)
	v.set_interior_value(v1)
	p.set_interior_value(p1)
	rho.set_interior_value(rho1)
	E.set_interior_value(E1)
	e.set_interior_value(e1)

	# Boundary Conditions:
	u.set_boundary_value( np.ones(ny)*( (y >= 0.5)*u0[1] + (y < 0.5)*u0[2] ), 'Dirichlet', 'Left')
	u.set_boundary_value( np.ones(ny)*( (y >= 0.5)*u0[0] + (y < 0.5)*u0[3] ), 'Dirichlet', 'Right')
	u.set_boundary_value( np.ones(nx)*( (x >= 0.5)*u0[3] + (x < 0.5)*u0[2] ), 'Dirichlet' 'Bottom')
	u.set_boundary_value( np.ones(nx)*( (x >= 0.5)*u0[0] + (x < 0.5)*u0[1] ), 'Dirichlet' 'Top')

	v.set_boundary_value( np.ones(ny)*( (y >= 0.5)*v0[1] + (y < 0.5)*v0[2] ), 'Dirichlet', 'Left')
	v.set_boundary_value( np.ones(ny)*( (y >= 0.5)*v0[0] + (y < 0.5)*v0[3] ), 'Dirichlet', 'Right')
	v.set_boundary_value( np.ones(nx)*( (x >= 0.5)*v0[3] + (x < 0.5)*v0[2] ), 'Dirichlet' 'Bottom')
	v.set_boundary_value( np.ones(nx)*( (x >= 0.5)*v0[0] + (x < 0.5)*v0[1] ), 'Dirichlet' 'Top')

	p.set_boundary_value( np.ones(ny)*( (y >= 0.5)*p0[1] + (y < 0.5)*p0[2] ), 'Dirichlet', 'Left')
	p.set_boundary_value( np.ones(ny)*( (y >= 0.5)*p0[0] + (y < 0.5)*p0[3] ), 'Dirichlet', 'Right')
	p.set_boundary_value( np.ones(nx)*( (x >= 0.5)*p0[3] + (x < 0.5)*p0[2] ), 'Dirichlet' 'Bottom')
	p.set_boundary_value( np.ones(nx)*( (x >= 0.5)*p0[0] + (x < 0.5)*p0[1] ), 'Dirichlet' 'Top')

	rho.set_boundary_value( np.ones(ny)*( (y >= 0.5)*rho0[1] + (y < 0.5)*rho0[2] ), 'Dirichlet', 'Left')
	rho.set_boundary_value( np.ones(ny)*( (y >= 0.5)*rho0[0] + (y < 0.5)*rho0[3] ), 'Dirichlet', 'Right')
	rho.set_boundary_value( np.ones(nx)*( (x >= 0.5)*rho0[3] + (x < 0.5)*rho0[2] ), 'Dirichlet' 'Bottom')
	rho.set_boundary_value( np.ones(nx)*( (x >= 0.5)*rho0[0] + (x < 0.5)*rho0[1] ), 'Dirichlet' 'Top')

	E.set_boundary_value( np.ones(ny)*( (y >= 0.5)*E0[1] + (y < 0.5)*E0[2] ) , 'Dirichlet', 'Left')
	E.set_boundary_value( np.ones(ny)*( (y >= 0.5)*E0[0] + (y < 0.5)*E0[3] ) , 'Dirichlet', 'Right')
	E.set_boundary_value( np.ones(nx)*( (x >= 0.5)*E0[3] + (x < 0.5)*E0[2] ) , 'Dirichlet' 'Bottom')
	E.set_boundary_value( np.ones(nx)*( (x >= 0.5)*E0[0] + (x < 0.5)*E0[1] ) , 'Dirichlet' 'Top')

	e.set_boundary_value( np.ones(ny)*( (y >= 0.5)*e0[1] + (y < 0.5)*e0[2] ) , 'Dirichlet', 'Left')
	e.set_boundary_value( np.ones(ny)*( (y >= 0.5)*e0[0] + (y < 0.5)*e0[3] ) , 'Dirichlet', 'Right')
	e.set_boundary_value( np.ones(nx)*( (x >= 0.5)*e0[3] + (x < 0.5)*e0[2] ) , 'Dirichlet' 'Bottom')
	e.set_boundary_value( np.ones(nx)*( (x >= 0.5)*e0[0] + (x < 0.5)*e0[1] ) , 'Dirichlet' 'Top')

	return [u,v,rho,p,E,e]

def Init_Bound_Variables2(u,v,rho,p,E,e,Mesh_Info,Field_Constants):
	'''
	This function initialises the variables and sets the boundary condition
	'''
	x  = Mesh_Info['x']
	X  = Mesh_Info['X']

	y  = Mesh_Info['y']
	Y  = Mesh_Info['Y']

	nx = Mesh_Info['nx']
	ny = Mesh_Info['ny']

	gamma = Field_Constants['gamma']

	# Inital Conditions
	p0   = np.array([ 1.0 , 0.1])
	rho0 = np.array([ 1.0 , 0.125])
	u0   = np.array([ 0.0 , 0.0])
	v0   = np.array([ 0.0 , 0.0])
	E0   = p0/(gamma-1) + (rho0/2.0) * (u0**2.0 + v0**2.0)	
	e0   = E0/rho0 -  (u0**2.0 + v0**2.0)

	p1   = p0[0]*(X<=0) + p0[1]*(X>0)
 	rho1 = rho0[0]*(X<=0) + rho0[1]*(X>0)
 	u1   = X*0.0
	v1   = X*0.0
	E1   = p1/(gamma-1) + (rho1/2.0) * (u1**2.0 + v1**2.0)
	e1   = E1/rho1 -  (u1**2.0 + v1**2.0)

	u.set_interior_value(u1)
	v.set_interior_value(v1)
	p.set_interior_value(p1)
	rho.set_interior_value(rho1)
	E.set_interior_value(E1)
	e.set_interior_value(e1)

	# Boundary Conditions:
	u.set_boundary_value( np.ones(ny)*( 0.0 ), 'Dirichlet', 'Left')
	u.set_boundary_value( np.ones(ny)*( 0.0 ), 'Dirichlet', 'Right')
	u.set_boundary_value( np.ones(nx)*( 0.0 ), 'ZeroGradient', 'Bottom')
	u.set_boundary_value( np.ones(nx)*( 0.0 ), 'ZeroGradient', 'Top')

	v.set_boundary_value( np.ones(ny)*( 0.0 ), 'Dirichlet', 'Left')
	v.set_boundary_value( np.ones(ny)*( 0.0 ), 'Dirichlet', 'Right')
	v.set_boundary_value( np.ones(nx)*( 0.0 ), 'ZeroGradient', 'Bottom')
	v.set_boundary_value( np.ones(nx)*( 0.0 ), 'ZeroGradient', 'Top')

	p.set_boundary_value( np.ones(ny)*( p0[0] ), 'Dirichlet', 'Left')
	p.set_boundary_value( np.ones(ny)*( p0[1] ), 'Dirichlet', 'Right')
	p.set_boundary_value( np.ones(nx)*( 0.0 ), 'ZeroGradient', 'Bottom')
	p.set_boundary_value( np.ones(nx)*( 0.0 ), 'ZeroGradient','Top')

	rho.set_boundary_value( np.ones(ny)*( rho0[0] ), 'Dirichlet', 'Left')
	rho.set_boundary_value( np.ones(ny)*( rho0[1] ), 'Dirichlet', 'Right')
	rho.set_boundary_value( np.ones(nx)*( 0.0 ), 'ZeroGradient', 'Bottom')
	rho.set_boundary_value( np.ones(nx)*( 0.0 ), 'ZeroGradient', 'Top')

	E.set_boundary_value( np.ones(ny)*( E0[0] ) , 'Dirichlet', 'Left')
	E.set_boundary_value( np.ones(ny)*( E0[-1] ) , 'Dirichlet', 'Right')
	E.set_boundary_value( np.ones(nx)*( 0.0 )   , 'ZeroGradient', 'Bottom')
	E.set_boundary_value( np.ones(nx)*( 0.0 )   , 'ZeroGradient', 'Top')

	e.set_boundary_value( np.ones(ny)*( e0[0] ) , 'Dirichlet', 'Left')
	e.set_boundary_value( np.ones(ny)*( e0[-1] ) , 'Dirichlet', 'Right')
	e.set_boundary_value( np.ones(nx)*( 0.0 )   , 'ZeroGradient', 'Bottom')
	e.set_boundary_value( np.ones(nx)*( 0.0 )   , 'ZeroGradient', 'Top')

	return [u,v,rho,p,E,e]
