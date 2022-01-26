import random
import math
import numpy as np
from resource import ind322

def multiple(N_parts, behavior):
####################################################################################################################
# - Receive number of parts N_parts for the scenario of interest and the behavior named as 'behavior';
# - Return True if the inequality is violated and False otherwise;
####################################################################################################################

	
	#-----------BOUND PARAMETERS----------------------------------------------------------
	
	atol=10**-5
	bound = 1.0+10**-8


	#-----------COMPUTING INEQUALITY------------------------------------------------------------
	
	(PI, PII) = Computing_PIPII(N_parts, behavior)

	EI = 2*PI - 1
	EII = 2*PII - 1
	Ineq = ((EI)**2 + (EII)**2)

	#print(Ineq)

	#-------------VERIFYING INEQUALITY VIOLATION-------------------------------------------------------

	if(Ineq>=bound):

		return True

	else:
		return False


def Computing_PIPII(N_parts, p):

####################################################################################################################
# - Receive number of parts N_parts for the scenario and the behavior p;
# - Return the inequality parameters P_I and P_II computed;
####################################################################################################################


	N = (2*N_parts) # Behavior's size
	P_I = 0 
	P_II = 0

	for i in range(0,2**N):

	#---Computing the bit-string I with (out/in)puts for p(a...|x...), respective to 'i'
		k = np.binary_repr(i)
		I = str(0)*(N-len(k))
		I = I+k

	#---Sum module 2 of all outputs a_i
		String_a = I[:int(len(I)/2)]
		String_a = np.array(list(map(int, list(String_a))))
		Sum_String_a = sum(String_a)%2

	#---Sum module 2 of all inputs x_i less the last
		String_x = I[int(len(I)/2):]
		String_x = np.array(list(map(int, list(String_x))))
		String_xm1 = np.delete(String_x, len(String_x)-1) 
		Sum_String_xm1 = sum(String_xm1)%2

		#print(I[len(I)-1], String_a, String_xm1)
		#if(int(I[len(I)-1]) == 0):
		#	print(Sum_String_a, Sum_String_xm1)
	#---Computing success probability for the first and second boolean functions respectively	
		if(Sum_String_a == 0 and int(I[len(I)-1]) == 0 ):

			P_I = P_I + p[i]

		if(Sum_String_a == Sum_String_xm1 and int(I[len(I)-1]) == 1 ):
			P_II = P_II + p[i]
	

	#---Re-normalizing
	P_I = (1/2**(N_parts-1))*P_I
	P_II = (1/2**(N_parts-1))*P_II

	return P_I, P_II



def Reordering_behavior(in_p):
####################################################################################################################
# - Receive a strictly 322 behavior vector ordered as ind322() function;
# - Return the behavior with re-ordered indexes with binary ordering indexes;
####################################################################################################################


	N = len(in_p)
	out_p = np.zeros (N, dtype='float')

	#(in)output bit-string size (a...x...)
	n = int(math.log(N, 2)) 

	for i in range(0,N):

	#---Computing the bit-string I with (out/in)puts for p(a...|x...), respective to 'i'
		k = np.binary_repr(i)
		I = str(0)*(n-len(k))
		I = I+k

	#---Translating indexes behavior
		out_p[i] = in_p[ind322(int(I[0]),int(I[1]),int(I[2]),int(I[3]),int(I[4]),int(I[5]))]
	
	return out_p


def Distributed_behavior(behavior_array):

####################################################################################################################
#  Receive a list of multipartite behaviors and return the distributed extended behavior
####################################################################################################################


#---Computing main parameters
	Behavior_sizes = [len(behavior_array[i]) for i in range(0,len(behavior_array))] # Size of each behavior
	Partitions = [int(math.log(Behavior_sizes[i], 2)/2) for i in range(0,len(Behavior_sizes))] # Computing number of parts for every partition

	N = 2*sum(Partitions)# Number of variables (inputs plus outputs)

	p = np.zeros(2**N, dtype='float') # Extended behavior initiation


#---Start to compute each component product
	for i in range(0,2**N):

	#---Computing the bit-string I with (out/in)puts for p(a...|x...), respective to 'i'
		k = np.binary_repr(i)
		I = str(0)*(N-len(k))
		I = I+k

		Input_string = I[:int(len(I)/2)]
		Output_string = I[int(len(I)/2):]


	#---Highlight in and output bit-string for each partition, compute the respective decimal index and perform the components  product 
		Component_product = 1
		del_left = 0
		del_right = 0 

		for j in range(0,len(Partitions)):

			del_right = del_right + Partitions[j]

			index = int(Input_string[del_left:del_right]+Output_string[del_left:del_right], 2)

			del_left = del_right

			Component_product = Component_product*behavior_array[j][index]


		p[i] = Component_product

	return p
