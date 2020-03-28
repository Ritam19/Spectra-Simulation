# -*- coding: utf-8 -*-
##################################################################################################################################################################################################
# Input:																							 #
#																								 #
# To run the program, the following files are needed:																		 #
#	1. spectra_simulation.py : file containing the code																	 #
#	2. inp.py : file containing the parameters and their values																 #
#	3. energies_forces.dat : file containing the excitation energies and oscillation strengths												 #
#          	                (to be modified by the user)																	 #
#	4. gamma.dat : file containing values of gamma, which will be used in some cases													 #
#                               (the values should be modified by the user)															 #
#	5. delta.dat : file containing values of delta, which will be used in some cases													 #
#                               (the values should be modified by the user)															 #
#  																								 #
#																								 #
# All the files are prepared and given, only the values should be modified															 #
#																								 #
# The file called “inp.dat” contains the parameters (gamma, delta, n, Emin, Emax, dE).														 #
#  These parameters can be modified according to the case as presented below.															 #
#                          																					 #
#																								 #
# Each parameter has a default value as shown:																			 #
#    																								 #
# Parameters	Allowed values						Unit		Default values												 #
# Gamma	        >0 – use this value as Γi for all bands																		 #
#               0 – read Γi from a file named gamma.dat 	         eV		    0.3												 #
# delta		≥0 – use this value as δi for all bands																		 #
#		0 – read δ from a file named delta.dat			 eV		    0.0													 #
# n		>0 – refractive index n				          -		    1.0													 #
# Emin		>0 – minimum value of E					 eV		    2.0													 #
# Emax		>0 – maximum value of E					 eV		    10.0												 #
# dE		>0 – interval size of E					 eV		    0.1													 #
#																								 #	
# Where the parameters are explained above.																			 #
#																								 #
# Two conditions are applied for all parameter:																			 #
#	If the parameter exists in the inp.py file, its value is used for all the excitation energies.												 #
#	If the parameter does not exist (deleted) in the file, the default value is used.													 #
#																								 #
#																								 #
# For gamma:																							 #
#	If the parameter exists and the given value > 0, this value is used for all the bands													 #
#	If the parameter exists and the given value = 0, gamma will be read from a “gamma.dat” file where it should contain a list of values corresponding to each excitation energy.		 #
#																								 #
#																								 #
#																								 #
#																								 #
# For delta:																							 #
#	If the parameter exists and the given value > or = 0, its value is used for all the bands												 #
#	If the parameter exists and the given value = -1, delta be read from a “delta.dat” file where it should contain a list of values corresponding to each excitation energy.		 #
#																								 #
#  2. The file called “energies_forces.dat” contains two columns, Ei ( excitation energies) and fi (oscillation strengths). 									 #
#                                     The values should be modified by the user.														 #
#																								 #
#																								 #
#  3. To be able to run the code, numpy library is needed.																	 #
#																								 #
#																								 #
#  4. To plot σ(E) x E:																						 #
#	Using matplotlib, install the library, uncomment the plotting part (the last four lines in the code) and the first line in the code							 #
#	Or simply using a graphical program as gnuplot or Excel, leave the commenting part.													 #
#																								 #
#																								 #
#  5. Run the program by: python spectra_simulation.py																		 #
#																								 #
#																								 #
#																								 #
#																								 #
# Output:																							 #
#																								 #
#																								 #
# 1. After execution , a file called “Output.dat” is created containing four columns:														 #
#                    ( E(eV)	λ(nm)	       σ(Å2.mol-1)	       ε(M-1cm-1) )														 #
#																								 #
#																								 #
#  Where:																							 #
#  λ is the wavenumber corresponding to (1240/E) 																		 #
#  ε is the molar absorbance or extinction coefficient (ε = σ/3.82353x10-5) 															 #
#																								 #
#																								 #
# 2. After every execution of the program, the Output.dat file in overwritten. It is recommended to rename the output file for each execution							 #
#																								 #
##################################################################################################################################################################################################




#import matplotlib.pyplot as plt
import sys
import numpy as np
from inp import *


Ei=np.genfromtxt('energies_forces.dat',skip_header=1,unpack=True,usecols=range(1))
fi=np.genfromtxt('energies_forces.dat',skip_header=1,unpack=True,usecols=range(1))
gamma_default= 0.3*np.ones(len(Ei))
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
        Ei,fi = np.genfromtxt('energies_forces.dat',skip_header=1,unpack=True,usecols=range(2))  
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

#plt.plot(E,sig)
#plt.xlabel('E')
#plt.ylabel('σ(E)')
#plt.show()
