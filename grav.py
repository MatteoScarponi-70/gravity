# 2021 Matteo Scarponi
# Forward gravity solution for a n-sided 2D polygon of given density
# A new class DensePoly is defined:
#	- inherits from Polygon class for handling geometries
#	- density property
#	- grav field calculation at given locations

from shapely.geometry import Polygon

class DensePoly(Polygon):

	def __init__(self, shell, density, holes=None):
		super().__init__(shell,holes)
		self.density=density
		if not self.is_valid:
			raise ValueError('This is not a valid 2D geometry')

	@property
	def density(self):
		print('Get density...')
		return self._density
	
	@density.setter
	def density(self,value):
		print('Set density...')
		if value<0:
			raise ValueError('Minimum density allowed: 0 kg/m3')
		self._density=value
	pass

# Example 

def main():
	import numpy as np
	import matplotlib.pyplot as plt

	R=5
	angles=np.arange(0,2*np.pi,0.01)
	x=np.cos(angles)
	y=np.sin(angles)

	poly=DensePoly(
		shell=list(zip(x,y)),
		density=200)

if __name__=='__main__':
	main()
