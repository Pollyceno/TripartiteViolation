import numpy as np
from resource import Index322


def uffink(n, p45, pd,  pw):
##################################################################################################################
# - Receives parameters' precision n, the three boxes that define the slice of interest p45, pd,  pw
# - Returns the parameters greatest gamma, for each epsilon, of the region pabcDxyz = gamma*p45 +epsilon*pd 
#   + (1-gamma-epsilon)*pw, for which the multiple copies inequality (FI)**2 + (FII)**2 <= 16 is not violated;
##################################################################################################################

	EM=np.arange(0., (1+1/n), 1/n)
	GM= 1 - EM
	
#-----------BISECTION PARAMETERS---------------------------------------------------------------
	iterac = 0
	gammaU = 1 
	gammaB = 0
	gamma = (gammaU - gammaB)/2
	atol=10**-5
	bound = 16.0+10**-8
#-----------VARING PARAMETERS------------------------------------------------------------------
	ç = n+1
	for i in range(0,ç):
		epsilon=(1/n)*i
		EM[i]=epsilon
		max_it = n
		aux = 0
		gammaU = 1-epsilon
		gamma = gammaU
		while(aux < max_it ):
			iterac = iterac +1
		
		#-----------SLICE----------------------------------------------------------------------
			pabcDxyz = gamma*p45 +epsilon*pd + (1-gamma-epsilon)*pw
			C001 = 0
			C010 = 0
			C100 = 0
			C111 = 0
			C110 = 0
			C101 = 0
			C011 = 0
			C000 = 0

			for a in range(0,2):
				for b in range(0,2):
					for c in range(0,2):
						
						C001 = C001 + ((-1)**(a+b+c))*pabcDxyz[Index322((a,b,c,0,0,1))]
						C010 = C010 + ((-1)**(a+b+c))*pabcDxyz[Index322((a,b,c,0,1,0))]
						C100 = C100 + ((-1)**(a+b+c))*pabcDxyz[Index322((a,b,c,1,0,0))]
						C111 = C111 + ((-1)**(a+b+c))*pabcDxyz[Index322((a,b,c,1,1,1))]

						C110 = C110 + ((-1)**(a+b+c))*pabcDxyz[Index322((a,b,c,1,1,0))]
						C101 = C101 + ((-1)**(a+b+c))*pabcDxyz[Index322((a,b,c,1,0,1))]
						C011 = C011 + ((-1)**(a+b+c))*pabcDxyz[Index322((a,b,c,0,1,1))]
						C000 = C000 + ((-1)**(a+b+c))*pabcDxyz[Index322((a,b,c,0,0,0))]
			#--0
			FI = C001 + C010 + C100 - C111
			FII = C110 + C101 + C011 - C000


		#-------------TESTING INEQUALITY-------------------------------------------------------
			desig = ((FI)**2 + (FII)**2)

			if(desig>bound):
				gammaU = gamma
				if(abs(desig - bound) > atol):#Up to the precision atol, devide the current gamma's range in half
					gamma = gammaU - (gammaU - gammaB)/2
				else:
					break
			else:
				gammaB = gamma
				if(abs(desig - bound) > atol):#Up to the precision atol, devide the current gamma's range in half
					gamma = gammaB + (gammaU - gammaB)/2
				else:
					break
			aux = aux +1

		aux = 0

		gammaB = 0.0 #reset
		GM[i] = gamma # Save the best gamma
	
	return (EM, GM)