"""PyWENO smooth reconstruction example."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
#from WENO import *

Col = ['C0','C1','C2','C3','C4','C5','C6','C7','C8','C9','b','g','m', 'r','y','#92C6FF','#E24A33', '#001C7F', '.15', '0.40','#30a2da', '#4C72B0', '#8EBA42', '#6d904f', '#bcbcbc', '#7600A1', '#D55E00', '#006374', '#B8860B', '#7A68A6', '#03ED3A', '0.70', '0.50', '#D65F5F', '#6ACC65', '#cbcbcb', '#55A868','#988ED5', '0.75','firebrick','#fc4f30','#E5E5E5', 'blue', '#fdb462', 'gray','#467821', '#feffb3', '#8b8b8b','#fa8174', '#64B5CD','#97F0AA', '#8C0900', 'm', '.8', '#8dd3c7', '#bfbbd9', 'green', '#8A2BE2', '0.60', '#E8000B', '#8172B2', '#FFC400', 'darkgoldenrod', '#CCB974', '#00FFCC','y', '#017517', 'c', '#D0BBFF', '#009E73', '#81b1d2', '#003FFF', '#b3de69', '0.00', 'purple', '#afeeee', '#77BEDB', '#CC79A7', '#348ABD', '#00D7FF', '#B0E0E6', '#EEEEEE','#ccebc4','b', '#555555', 'k', '#bc82bd', 'r', '#C4AD66', '#F0E442', '#777777', '#ffed6f', '#B47CC7', '#C44E52', '#eeeeee', '#FBC15E', '#56B4E9', '#4878CF', '#FFFEA3', '#f0f0f0', '#FFB5B8', '#e5ae38', '#FF9F9A', '#EAEAF2', '0.5', 'black', 'g', '#A60628','#0072B2','red']

Mar = ["x","s","o","v","^","<",">","d","2","3","4","1","+","p","P","*","h","H","x","D","d","|","_",""]

Lin = ['-','--',':','-.', (0, (3, 1, 1, 1, 1, 1)), (0, (5, 10)), (0, (3, 10, 1, 10)), (0, (3, 10, 1, 10, 1, 10))]

def plot_fig(x,u,t,k,Name,limits,Fig_Name,ylabel_name='u(x)',marker=True,x_opt=[],y_opt=[],name_opt=None,Leg_loc = 'upper right'):
	L = ['-',':','-.','--']
	#fig, axes = plt.subplots(3,2,figsize=(8, 6)) #, figsize=(8, 4))
	fig, axes = plt.subplots(figsize=(8, 6)) #, figsize=(8, 4))
	ax = axes

	if marker==True:
		for i in range(len(u)):
			ax.plot(x,u[i],label=Name[i],Marker=Mar[i],linestyle=L[i%4])
		if name_opt != None:
			i +=1
			ax.plot(x_opt,y_opt,label=name_opt,Marker=Mar[i],linestyle=L[i%4])
	else:
		for i in range(len(u)):
			ax.plot(x,u[i],label=Name[i])
		if name_opt != None:
			i +=1
			ax.plot(x_opt,y_opt,label=name_opt,linestyle=L[i%4])

	ax.set_ylabel(ylabel_name)
	ax.set_xlabel('x')

	#ax.set_xticks([])
	#ax.set_xticklabels([],fontweight=0)
	ax.legend(fontsize=10,loc=Leg_loc)
	ax.set_xlim([limits[0],limits[1]])
	ax.set_ylim([limits[2],limits[3]])


	fig.tight_layout()
	fig.savefig(Fig_Name+'-%d.png' % t )
	print 'figure saved in ',Fig_Name+'-%d.png' % t 
	plt.close()



def plot_contour3D(X,Y,u,t,Name,limits):

	from mpl_toolkits.mplot3d import Axes3D

	fig = plt.figure(figsize=(11, 7), dpi=100)
	ax = fig.gca(projection='3d')
	ax.plot_surface(X, Y, u[:], cmap=cm.viridis, rstride=1, cstride=1)
	ax.set_xlabel('$x$')
	ax.set_ylabel('$y$')
	ax.set_zlim(limits[0], limits[1])
	plt.savefig(Name)
	plt.close()



def plot_contour(X,Y,u,t,Name,limits):

	fig, ax = plt.subplots(1) 
	levels = np.linspace(limits[0], limits[1],21)

	ax.contourf(X, Y, u,levels,cmap = cm.viridis,extend='both')
	ax.set_xlabel('$x$')
	ax.set_ylabel('$y$')
	plt.savefig(Name)
	plt.close()
