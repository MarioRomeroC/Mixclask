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
        self.__target_wl = target_wavelength
        self.__max_it    = max_iterations
        self.tolerance   = tolerance
        
        self.n_iterations = 0
        
        self.__prev_values = [None for i in range(0,len(self.__sedFiles))]
        self.__prev_result = None
        
    def __compute_values(self):
        result = np.empty(len(self.__sedFiles))
        for ii in range(0,len(self.__sedFiles)):
            file = self.__sedFiles[ii]
            all_data = unk.readColumn(file,[0,1])
            wavelengths  = all_data[:,0]
            FourPi_nuJnu = all_data[:,1]
            
            #For now, I get the value at the target wavelength to check convergence.
            result[ii] = np.interp(self.__target_wl,wavelengths,FourPi_nuJnu)
            
        return result
    
    def stop(self):
        # Some updates first
        curr_values = self.__compute_values()
        self.n_iterations += 1
        
        #Now check things
        if self.n_iterations > self.__max_it:
            #Stop if maximum number of iteration has been reached
            print("Stopped because max iterations reached")
            return True
        elif all(elem != None for elem in self.__prev_values):
            curr_result = unk.RMSD(self.__prev_values,curr_values) #root mean squared deviation
            error = None
            if self.__prev_result != None:
                diff  = abs(curr_result-self.__prev_result)
                error = (diff/self.__prev_result) if diff != 0.0 else 0.0 #To avoid a (0-0)/0 = nan and so never converges
            else:
                error = self.tolerance*1000.0 #1000.0 because potato, I only want to force next 'else' sentence
            #print(error) #Debug
            if error < self.tolerance:
                print("Stopped because convergence has been achieved")
                return True
            else:
                self.__prev_result = curr_result
        
        self.__prev_values = curr_values
        
        return False
    
    def iteration(self,time_elapsed):
        print("iteration "+str(self.n_iterations)+" ("+str(time_elapsed)+" s) :")
        return self.n_iterations