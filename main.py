import numpy as np
import os
import time

from resource import resource
#from npa322 import NPA322
from multiple import *
from group import *


###########################################################################################################################
			# BOXES
###########################################################################################################################


Class = []
t = 0

for N in range(1,47):
	
	#pabcDxyz = recurso2() #Generates the N class representative
	pabcDxyz = resource(N)

	resources = All_permutations(pabcDxyz) # Generates all class's permutation
	print('Class: ',N, 'Size: ', len(resources))
	box = 0
	L=False
	
	for p in resources:
		#print('Box: ', box)

		Dist_p = Distributed_behavior(np.array([Reordering_behavior(p), Reordering_behavior(p), Reordering_behavior(p)]))

		T = time.time()
		L = multiple(9, Dist_p) # Test violation
		#L = multiple(3, Reordering_behavior(p))
		
		
		if(L==True):
			print('Box ',box, L)
			
			Class = np.loadtxt('Violating_classes.txt', dtype='int')
			
			Class = np.sort(np.unique(np.append(Class, N)))

			np.savetxt('Violating_classes.txt', Class, fmt="%s")

			break
		
		print("Spent time--- %s seconds ---" % (time.time() - T))

		box = box +1

		#break

	if(L==False):
		print('Does not')

	break
	


