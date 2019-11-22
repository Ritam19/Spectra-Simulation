# -*- coding: utf-8 -*-
#################################################################################################################################################################################################
#Simple Spectra simulation program                                                                                                                                                              #
#(released 25.10.2019)                                                                                                                                                                          #
#                                                                                                                                                                                               #
#How  to run it?                                 																		#
#																								#
#Input:																						                #
#All the file are given as an example, only the values should be modified								    			    				#
#1-An input file called “inp.py” is prepared containing all the parameters:															#
#		(gamma, delta, n, Emin, Emax, dE)																		#
#																								#
#Each parameter has a default value as shown:																		        #
#    																								#
#Parameters	Unit	Default value																				#
#gamma		eV	0.4																					#
#delta		eV	0.0																					#
#n		-	1.0																					#
#Emin		eV	2.0																					#
#Emax		eV	10.0																					#
#dE		eV	0.1																					#
#																								#
#																								#
#Where  gamma is the width of a band i																				#
#       delta is the energy shift between the vertical excitation and the band maximum														#
#	n is the refractive index																				#
#	Emin is the minimum value of E																				#
#       Emax is the maximum value of E																				#
#       dE is the interval size of E																				#
#																								#
#Two conditions are applied for all parameter:																			#
#a)	If the parameter exists in the inp.py file, its value is used for all the excitation energies.												#
#b)	If the parameter does not exist (deleted) in the file, the default value is used.										    			#
#																								#
#																								#
#For gamma:																							#
#a)	If the parameter exists and the given value > 0, this value is used for all the bands											    		#
#b)	If the parameter exists and the given value = 0, gamma should be read from a “gamma.dat” file where it contains a list of values coresponding to each excitation energy.    	    	#
#																								#
#																								#
#For delta:																							#
#a)	If the parameter exists and the given  value > or = 0, this value is used for all the bands											    	#
#b)	If the parameter exists and the given value = -1, delta should be read from a “delta.dat” file where it contains a list of values coresponding to each excitation energy.   	    	#
#																								#
#																								#
#2-energies_forces.py file should be modified, containing two list (Ei and fi) 	                                                                                                		#
#																								#
#																								#
#3-Run the program by: python spectra_simulation.py																	    	#
#																								#
#																								#
#																								#
#																								#	
#																								#
#Output:																							#
# 																								#
#1-After execution , an “Output.dat” file is created containing four columns:															#
#                    ( E(eV)	   λ(nm)	       σ(Å2.mol-1)        ε(M-1cm-1) )														#
#																								#
#Where λ is the wavenumber corresponding to (1240/E) , and ε is the molar absorbance or extinction coefficient (ε = σ/3.82353x10-5) 								#
#																								#
#																								#
#2-To plot (E) x σ(E): 																					        #
#  a)install "matplotlib" library, and uncomment the last four lines in the code														#
#  b)use a graphical porgram as gnuplot or Excel																		#
#																								#
#################################################################################################################################################################################################



import matplotlib.pyplot as plt
import sys
import numpy as np
from inp import *

Ei=np.genfromtxt('energies_forces.py',skip_header=1,unpack=True,usecols=range(1))
fi=np.genfromtxt('energies_forces.py',skip_header=1,unpack=True,usecols=range(1))
gamma_default= 0.4*np.ones(len(Ei))
delta_default=0.0*np.ones(len(Ei))
nref_default=1.0
Emin_default=2.0
Emax_default=10.0
step_default=0.1

#-----------------------------#
# Calculation of the energies #
#-----------------------------#

E = np.arange(Emin, Emax, step)



#-------------------------------#
# Calculation of the wavenumber #
#-------------------------------#

lam=np.divide(1240,E)



#--------------------------------#
# Function to calculate spectrum #
#--------------------------------#

def calculate_spectra(Ei,fi,E,gamma,delta,nref):
        Ei,fi = np.genfromtxt('energies_forces.py',skip_header=1,unpack=True,usecols=range(2))  
        iE=len(E)
        sig = np.zeros(iE)
        for j in range(iE):
                res=0
                for i in range(0,len(Ei)):      
                        a=np.multiply(0.619,nref)
                        z=(-((E[j]-Ei[i]+delta[i])**2)/(gamma[i]**2))
                        res=res+((fi[i]/gamma[i])*np.exp(z)) 
                sig[j]=(a*res)
        return(sig)


#--------------#
# MAIN PROGRAM #
#--------------#
def calc():
        gamma_c = 0
        delta_c = 0
        nref_c = 0
        Emin_c = 0
        Emax_c = 0
        step_c = 0
        infile = open('inp.py', 'r')
        for line in infile:
                variable, value = line.split('=')
                variable = variable.strip()  
                value=float(value)
                if gamma_c == 0 and variable == 'gamma':
                        if value == 0: 
                            gamma = np.genfromtxt('gamma.dat',skip_header=1,usecols=range(1))
                        elif value < 0:
                          print('Invalid value')
                          sys.exit()
                        elif value > 0:
                          gamma = value*np.ones(len(Ei))
                        gamma_c = 1
                elif gamma_c == 1:
                        pass
                else:
                        gamma=gamma_default
                if delta_c == 0 and variable == 'delta':
                        if value == -1: 
                            delta=np.genfromtxt('delta.dat',skip_header=1,usecols=range(0))
                        elif value == 0:
                          print('Invalid value')
                          sys.exit()
                        elif value >= 0:
                          delta = value*np.ones(len(Ei))
                        delta_c = 1
                elif delta_c ==1:
                        pass
                else:   
                        delta=delta_default
                if nref_c == 0 and variable == 'nref':
                        nref = float(value)
                        nref_c = 1
                elif nref_c ==1:
                        pass
                else:
                        nref=nref_default
                if Emin_c == 0 and variable == 'Emin':
                        Emin = float(value)
                        Emin_c = 1
                elif Emin_c ==1:
                        pass
                else:
                        Emin=Emin_default
                if Emax_c == 0 and variable == 'Emax':
                        Emax = float(value)
                        Emax_c = 1
                elif Emax_c == 1:
                        pass
                else:
                        Emax=Emax_default
                if step_c == 0 and variable == 'step':
                        step = float(value)
                        step_c = 1
                elif step_c == 1:
                        pass    
                else:
                        step=step_default

        return gamma,delta,nref,Emin,Emax,step
        infile.close()
         
gamma,delta,nref,Emin,Emax,step = calc()

#---------------#
# Call function #
#---------------#

sig= calculate_spectra(Ei,fi,E,gamma,delta,nref)


#---------------------#
# Calculating epsilon #
#---------------------#

epsi = np.divide(sig,(3.82353*10**-5))


#---------#
# Outputs #
#---------#

energies=E
lamda=lam
sigma=sig
epsilon=epsi


#------------------------------#
# Saving the Outputs in a file #
#------------------------------#

Output=np.column_stack((energies,lamda,sigma,epsilon))
header= "E(eV)          λ(nm)                σ(Å2.mol-1)      ε(M-1cm-1)\n"
results=np.savetxt('Output.dat', Output, fmt='%.5f', delimiter= "\t\t", header=header)




#-----------#
# plotting #
#-----------#

plt.plot(E,sig)
plt.xlabel('E')
plt.ylabel(r'$\sigma$')
plt.show()
