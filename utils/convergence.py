#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class and methods needed in order to check convergence

Mario Romero        November 2021
"""

import numpy as np
import utils.unkeep as unk

class ConvergenceObject(object):
    #Not defined as function because I need to store parameters
    def __init__(self,SED_files,target_wavelength,max_iterations,tolerance):
        #Declare some parameters
        self.__sedFiles  = SED_files #like CloudyClass
        self.__max_it    = max_iterations
        self.tolerance   = tolerance
        
        self.n_iterations = 0
        
        if isinstance(target_wavelength,(list,tuple,np.ndarray)):
            self.__target_wl = target_wavelength #Can be a list
        else: #int or float
            self.__target_wl = [target_wavelength]
        
        #Store previous values:
        self.__prev_values = [ [ None for j in range(0,len(self.__sedFiles)) ] for i in range(0,len(self.__target_wl))]
        self.__prev_result = [ None for i in range(0,len(self.__target_wl)) ]
        #so:
        # self.__prev_values[wavelength][zone]
        # self.__prev_result[zone]
    
    def __readSEDfile(self,file):
        #SEDfiles have their data sorted from highest wavelength to shortest. Numpy wants the opposite order, thus np.flip solves this issue
        all_data = unk.readColumn(file,[0,1])
        wavelengths  = np.flip(all_data[:,0])
        FourPi_nuJnu = np.flip(all_data[:,1])
        
        return wavelengths, FourPi_nuJnu
    
    def __compute_values(self,desired_wl):
        if isinstance(desired_wl,(tuple,list,np.ndarray)):
            return self.__compute_integral(desired_wl)
        else: #you gave a float
            return self.__find_values(desired_wl)
    
    def __find_values(self,desired_wl):
        #desired_wl is a FLOAT
        result = np.empty(len(self.__sedFiles))
        for ii in range(0,len(self.__sedFiles)):
            #file = self.__sedFiles[ii]
            #all_data = unk.readColumn(file,[0,1])
            #wavelengths  = np.flip(all_data[:,0])
            #FourPi_nuJnu = np.flip(all_data[:,1])
            wavelengths, FourPi_nuJnu = self.__readSEDfile(self.__sedFiles[ii])
            
            #For now, I get the value at the target wavelength to check convergence.
            result[ii] = np.interp(desired_wl,wavelengths,FourPi_nuJnu)
            
        return result
    
    def __compute_integral(self,desired_wl):
        #desired_wl is a tuple/list -> (x[0],x[1])
        result = np.empty(len(self.__sedFiles))
        for ii in range(0,len(self.__sedFiles)):
            wavelengths, FourPi_nuJnu = self.__readSEDfile(self.__sedFiles[ii])
            result[ii] = unk.integrate(desired_wl,wavelengths,FourPi_nuJnu)
        
        return result
            
    
    def iteration(self,time_elapsed):
        print("iteration "+str(self.n_iterations)+" ("+str(time_elapsed)+" s) :")
        return self.n_iterations
    
    def stop(self):
        self.n_iterations += 1
        #First, check if we have reached max_iterations
        if self.n_iterations > self.__max_it:
            #Stop if maximum number of iteration has been reached
            print("Warning: Stopped because max iterations reached")
            return True
        
        #Now check, per target_wavelength, if error < tolerance
        convergence = []
        error = []
        for w in range(0,len(self.__target_wl)):
            curr_values_w = self.__compute_values(self.__target_wl[w]) #Here, we look for the value of 4pi*lambda*J_lambda at the wavelength
            if all(elem == None for elem in self.__prev_values[w]): 
                #We have not iterated before, so keep this result
                #This happens for the iteration 0 (gas is not included yet)
                convergence.append(False)
                error.append(None)
            elif self.__prev_result[w] == None:
                #We have iterated once, and because computing the error means to have current and previous result, this value has not been computed before
                #This happens for iteration 1
                self.__prev_result[w] = unk.RMSD(self.__prev_values[w],curr_values_w) #root mean squared deviation
                convergence.append(False)
                error.append(None)
            else:
                #We have iterated twice. We can now check convergence
                curr_result_w = unk.RMSD(self.__prev_values[w],curr_values_w)
                diff = abs(curr_result_w-self.__prev_result[w])
                error.append( (diff/self.__prev_result[w]) if diff != 0.0 else 0.0 )  #To avoid a (0-0)/0 = nan and so never converges
                if error[-1] < self.tolerance:
                    convergence.append(True)
                else:
                    convergence.append(False)
                
                self.__prev_result[w] = curr_result_w
                
            self.__prev_values[w] = curr_values_w
        
        
        print("Debug: Tolerance = "+str(self.tolerance)+" . Errors = "+str(error)+" . Wavelengths = "+str(self.__target_wl))
        if any(elem == False for elem in convergence):
            print("Convergence not reached, I will iterate again.")
            return False
        else:
            print("Convergence achieved, stopping...")
            return True
