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
        self.__prev_intensity = [ [ None for j in range(0,len(self.__sedFiles)) ] for i in range(0,len(self.__target_wl))]
        self.__prev_result = [ None for i in range(0,len(self.__target_wl)) ]
        #so:
        # self.__prev_intensity[wavelength][zone]
        # self.__prev_result[zone]
    
    def __readSEDfile(self,file):
        #SEDfiles have their data sorted from highest wavelength to shortest. Numpy wants the opposite order, thus np.flip solves this issue
        all_data = unk.readColumn(file,[0,1])
        wavelengths  = np.flip(all_data[:,0])
        FourPi_nuJnu = np.flip(all_data[:,1])
        
        return wavelengths, FourPi_nuJnu
    
    def __compute_intensity(self,desired_wl):
        if isinstance(desired_wl,(tuple,list,np.ndarray)):
            return self.__integrated_intensity(desired_wl)
        else: #you gave a float
            return self.__interpolated_intensity(desired_wl)
    
    def __interpolated_intensity(self,desired_wl):
        #desired_wl is a FLOAT
        result = np.empty(len(self.__sedFiles))
        for ii in range(0,len(self.__sedFiles)):
            wavelengths, FourPi_nuJnu = self.__readSEDfile(self.__sedFiles[ii])
            
            #Get the value at the target wavelength to check convergence.
            result[ii] = np.interp(desired_wl,wavelengths,FourPi_nuJnu)
            
        return result
    
    def __integrated_intensity(self,desired_wl_range):
        #desired_wl is a tuple/list -> (x[0],x[1])
        result = np.empty(len(self.__sedFiles))
        for ii in range(0,len(self.__sedFiles)):
            wavelengths, FourPi_nuJnu = self.__readSEDfile(self.__sedFiles[ii])
            result[ii] = unk.integrate(desired_wl_range,wavelengths,FourPi_nuJnu)
        
        return result
            
    
    def iteration(self,time_elapsed):
        print("iteration "+str(self.n_iterations)+" ("+str(time_elapsed)+" s) :")
        return self.n_iterations
    
    def stop(self):
        self.n_iterations += 1
        #First, check if we have reached max_iterations
        if self.n_iterations > self.__max_it:
            #Stop if maximum number of iteration has been reached
            print("Warning: Stopped because max iterations reached. Check if tolerance is too low for the montecarlo noise")
            return True
        
        #Now check, per target_wavelength, if error < tolerance
        convergence = []
        for w in range(0,len(self.__target_wl)):
            print("Debug: current wavelength = "+str(self.__target_wl[w]))
            curr_intensity = self.__compute_intensity(self.__target_wl[w])
            if all(elem == None for elem in self.__prev_intensity[w]): 
                #We have not iterated before, so keep this result
                #This happens for the iteration 0 (gas is not included yet)
                convergence.append(False)
            else:
                #We have iterated once. We can now check convergence
                intensity_difference = abs(curr_intensity - self.__prev_intensity[w])
                intensity_mean       = 0.5 * (curr_intensity + self.__prev_intensity[w])
                n_zones = len(intensity_difference) #Just for clarity
                if all(intensity_difference[zone] <= self.tolerance*intensity_mean[zone] for zone in range(0,n_zones)):
                    convergence.append(True)
                else:
                    convergence.append(False)
            
                #print("Debug: Differences = "+str(intensity_difference)+" . Means = "+str(intensity_mean))
            
            self.__prev_intensity[w] = curr_intensity
        
        
        #print("Debug: Tolerance = "+str(self.tolerance)+" . Convergence? : "+str(convergence))
        if any(elem == False for elem in convergence):
            print("Convergence not reached, I will iterate again.")
            return False
        else:
            print("Convergence achieved, stopping...")
            return True
