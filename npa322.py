import ncpol2sdpa as nc
from recurso import *
from NCP import *

def NPA322(n, p44, pd,  pw, level):
		
	E=np.zeros(n+1, dtype='float')
	G=np.zeros(n+1, dtype='float')

##################################################################################################
##################### PARAMETROS BISECCAO ########################################################

	#-----------PARAMETROS BISECCAO----------------------------------------------------------
	iterac = 0
	gammaU = 1 
	gammaB = 0
	gamma = (gammaU - gammaB)/2
	atol=0.001
	bound = 0.0
	#deltaU = 0.01
	#deltaB = 0.15
##################################################################################################
##################### VARIANDO PARAMETROS ########################################################

	ç = n+1
	#ç=n 
	for i in range(0,ç):
		epsilon=(1/n)*i
		E[i]=epsilon
		max_it = 15
		aux = 0
		gammaU = 1-epsilon
		gamma = gammaB
		#print(gamma, gammaU, epsilon)
		while(aux < max_it ):
			iterac = iterac +1

##################################################################################################
##################### RECURSO ####################################################################

			pabcDxyz = gamma*p44 + (1-gamma)*pw
			#gamma = 0.0
			#pabcDxyz = gamma*p44 + epsilon*pd + (1-gamma-epsilon)*pw


##################################################################################################
##################### VERIFICANDO  ###############################################################
			Q = quantumedge(pabcDxyz, level)


			if(Q == 'dual_infeas_cer'):
				print('U',epsilon, gamma, Q)
				gammaU = gamma
				if(abs(gammaU - gammaB) > atol):
					gamma = gammaU - (gammaU - gammaB)/2
				else:
					break
			else:
				print('B',epsilon, gamma, Q)
				gammaB = gamma
				if(abs(gammaU - gammaB) > atol):
					gamma = gammaB + (gammaU - gammaB)/2
				else:
					break
			aux = aux +1
	
		aux = 0
	
		gammaB = 0.0 #gamma - deltaB
		G[i] = gamma




	return (E, G)