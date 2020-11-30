import numpy as np

def E_total_energy(rho,u,v,p,gamma):
	return rho*0.5*(u**2.0 + v**2.0) + p/(gamma-1) # E = rho*(e + 0.5*(u^2 + v^2))

def p_idealgas(E,rho,u,v,gamma):
	return (gamma-1)*(E - rho*0.5*(u**2.0+v**2.0) )

def e_internal_energy(E,rho,u,v):
	return (E/rho) - 0.5*(u**2.0 + v**2.0)