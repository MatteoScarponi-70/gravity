# 2021 Matteo Scarponi
# Forward gravity solution for a n-sided 2D polygon of given density

# Implemented formulas are published by Won and Bevis 1987
#	Won, I.J. and Bevis, M., 1987. Computing the gravitational 
# 	and magnetic anomalies due to a polygon: Algorithms and 
#	Fortran subroutines. Geophysics, 52(2), pp.232-238.

# Here, class "DensePoly" is defined:

#	- inherits from shapely.geometry.Polygon class
#	- defines density as new class property
#	- computes gravity signal at given locations
#	- distances should be given in meters [m]
#	- density should be given in [kg/m3]

from shapely.geometry import Polygon
import numpy as np

class DensePolygon(Polygon):

	def __init__(self, shell, density, holes=None):
		super().__init__(shell,holes)
		self.density=density
		if not self.is_valid:
			raise ValueError('This is not a valid 2D geometry')

	@property
	def density(self):
		return self._density
	
	@density.setter
	def density(self,value):
		if value<0:
			raise ValueError('Minimum density allowed: 0 kg/m3')
		self._density=value

	def Grav(self,xstations,zstations):

		G=6.67259e-11 # Gravitational Constant m3/kg/s2

		xPolygon,zPolygon=self.exterior.xy
		xPolygon=np.array(xPolygon)
		zPolygon=np.array(zPolygon)

		zAnomaly=np.zeros(xstations.shape)
		xAnomaly=np.zeros(xstations.shape)

		for i in range(len(xstations)):

			print('Station '+str(i),end='\r',flush=True)

			x=xPolygon-xstations[i]
			z=zPolygon-zstations[i]

			z[z==0]=1e-42

			theta=np.arctan2(z,x)
			r=np.sqrt((x**2)+(z**2))

			X=np.zeros(len(xPolygon)-1)
			Z=np.zeros(len(xPolygon)-1)

			for j in range(len(xPolygon)-1):

				theta1=theta[j]
				theta2=theta[j+1]

				# Case 1

				if z[j]*z[j+1]<0:

					if (x[j]*z[j+1]) < (x[j+1]*z[j]) and z[j+1]>=0:
						theta1=theta1+2*np.pi

					if (x[j]*z[j+1]) > (x[j+1]*z[j]) and z[j]>=0:
						theta2=theta2+2*np.pi

					if (x[j]*z[j+1]) == (x[j+1]*z[j]):
						Z[j]=0
						X[j]=0
						continue

				# Case 2

				if (x[j]==0 and z[j]==0) or (x[j+1]==0 and z[j+1]==0):
					Z[j]=0
					X[j]=0
					continue

				# Case 3

				if (x[j]==x[j+1]):
					Z[j]=x[j]*np.log(r[j+1]/r[j])
					X[j]=-x[j]*(theta1-theta2)
					continue

				# General case

				A=(x[j+1]-x[j])*(x[j]*z[j+1]-x[j+1]*z[j])/((x[j+1]-x[j])**2+(z[j+1]-z[j])**2)
				B=(z[j+1]-z[j])/(x[j+1]-x[j])

				Z[j]=A*((theta1-theta2)+B*np.log(r[j+1]/r[j]))
				X[j]=A*(B*(theta2-theta1)+np.log(r[j+1]/r[j]))

				continue

			zAnomaly[i]=2*G*self.density*np.sum(Z)*1e5 # mGal
			xAnomaly[i]=2*G*self.density*np.sum(X)*1e5 # mGal

		return xAnomaly,zAnomaly

def main():
	import numpy as np
	import matplotlib.pyplot as plt

	# Exemple #

	R=1		# km
	dz=5	# km
	angles=np.arange(0,2*np.pi,0.01)	# Radians
	x=R*np.cos(angles)
	z=R*np.sin(angles)+dz

	x=x*1000 # m
	z=z*1000 # m

	poly=DensePolygon(
		shell=list(zip(x,z)),
		density=200)

	xstations=np.arange(-25,25,0.5)*1000		# m
	zstations=np.zeros(xstations.shape)*1000	# m

	xAnomaly,zAnomaly=poly.Grav(xstations,-zstations) # mGal

	# Plot results #

	f=plt.figure(1,figsize=[8,6])
	ax=plt.subplot(212)
	ax.fill(
		*poly.exterior.xy,
		color='#F2A128',
		edgecolor='k',
		label='Body')
	ax.scatter(
		xstations,
		zstations,
		50,
		facecolors='r',
		edgecolors='k',
		marker='v',
		lw=0.2,
		label='Stations')
	ax.grid(
		True,
		ls='--',
		lw=0.5,
		alpha=0.4)
	ax.set_xlim(-25000,25000)
	ax.set_ylim(8000,0)
	ax.set_aspect('equal')
	ax.set_xlabel('Dist [m]')
	ax.set_ylabel('Depth [m]')
	ax.legend(loc='best')
	ax=plt.subplot(211,sharex=ax)
	ax.scatter(
		xstations,
		zAnomaly,
		50,
		edgecolors='k',
		facecolors='r',
		lw=0.2,
		label='gz')
	ax.scatter(
		xstations,
		xAnomaly,
		50,
		edgecolors='k',
		facecolors='g',
		lw=0.2,
		label='gx')
	ax.grid(
		True,
		ls='--',
		lw=0.5,
		alpha=0.4)
	ax.set_ylabel('Gravity anomaly [mGal]')
	ax.legend(loc='best')
	plt.tight_layout()
	plt.show()
	plt.close()

if __name__=='__main__':
	main()