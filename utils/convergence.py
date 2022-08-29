#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Class and methods needed in order to check convergence.
Architecture has been changed to check convergence in several quantities in the future
(instead of only radiation field, although now it is the only option).

Mario Romero        August 2022
"""

import numpy as np
import utils.unkeep as unk

class ConvergenceObject(object):
    def __init__(self,SED_files,options_dict):
        #Store some internal parameters first
        self.n_iterations = 0
        self.__start_iteration = options_dict['Technical']['iteration_to_start_convergence']
        self.__max_iterations = options_dict['Technical']['n_iterations']
        self.__pointers = SED_files #These are the filenames to extract info from it
        self.__n_zones = len(self.__pointers)
        self.__criteria = options_dict['Convergence']['Criteria']
        #
        self.__initStorage(options_dict['Convergence'])

    def __initStorage(self,conv_dict):
        '''
        This method converts the semi-structured data from conv_dict into structured lists.
        The output list is self.__convResults and its structure is:
            x[ KEY ][ ITERATION ][ ZONE ][ NUMBER ]
        NUMBER is:
        0 : Result computed with 'computeData' function below
        1 : Mean of 0 for all x[fixed_key][][fixed_zone][0] between iteration 1 and current iteration.
        2 : Variance of 0 for all x[fixed_key][][fixed_zone][0] between iteration 1 and current iteration.
        3 : Median of 0 for all x[fixed_key][][fixed_zone][0] between iteration 1 and current iteration.
        4 : x84-x16 of 0 for all x[fixed_key][][fixed_zone][0] between iteration 1 and current iteration (x84 and x16 are percentiles 84 and 16).
        Iteration 0 is skipped because there is no gas in it
        '''
        n_quantities = len(conv_dict)-1
        self.__names = [] #For testing purposes
        self.__tolerances = []
        self.__wlRanges = []
        for key in conv_dict:
            if key == 'Criteria': continue
            self.__names.append(key)
            self.__tolerances.append(conv_dict[key]['tolerance'])
            self.__wlRanges.append(conv_dict[key]['wavelengthRange'])

        self.__convResults = np.nan*np.ones((n_quantities,self.__max_iterations+1,self.__n_zones,5)) #I init to NaN to indentify not filled data
        #max_iterations+1 because you still have to consider iteration0 (the one without ISM)
        #That's all!

    def __readFile(self, file):
        # Files from SEDfiles have their data sorted from highest wavelength to shortest.
        # Numpy wants the opposite order, thus np.flip solves this issue
        all_data = unk.readColumn(file, [0, 1])
        wavelengths = np.flip(all_data[:, 0])
        function = np.flip(all_data[:, 1])

        return wavelengths, function

    ### PUBLIC METHODS ###
    def iteration(self,time_elapsed):
        print("iteration "+str(self.n_iterations)+" ("+str(time_elapsed)+" s) :")
        return self.n_iterations

    def stop(self):
        '''
        This is the main method of this class, it returns True if convergence has been achieved and vice-versa.
        It is called per iteration (so it is FIXED in this method, to be increased by one at the end)
        It also contains several subroutines to compute data.
        '''
        #Compute data function
        def computeData(curr_key,curr_zone):
            #Read the correct files
            wavelengths, FourPi_nuJnu = self.__readFile(self.__pointers[curr_zone])
            if isinstance(self.__wlRanges[curr_key],(tuple,list,np.ndarray)): #Wavelength range given
                #Integrate data between the wavelength range
                return unk.integrate(self.__wlRanges[curr_key], wavelengths, FourPi_nuJnu)
            else: #Fixed wavelength given
                #So interpolate
                return np.interp(self.__wlRanges[curr_key],wavelengths,FourPi_nuJnu)

        #Store data
        for zone in range(0,self.__n_zones):
            for key in range(0,len(self.__names)):
                curr_result = computeData(key,zone)
                self.__convResults[key][self.n_iterations][zone][0] = curr_result
                self.__convResults[key][self.n_iterations][zone][1] = np.nanmean(self.__convResults[key,:,zone,0])
                self.__convResults[key][self.n_iterations][zone][2] = np.nanvar(self.__convResults[key, :, zone, 0])
                self.__convResults[key][self.n_iterations][zone][3] = np.nanmedian(self.__convResults[key, :, zone, 0])
                self.__convResults[key][self.n_iterations][zone][4] = abs( np.nanpercentile(self.__convResults[key, :, zone, 0],84)
                                                                           - np.nanpercentile(self.__convResults[key, :, zone, 0],16) )

        #Check convergence
        has_converged = None
        #To have a True as a result, you need that ALL ZONES PER ALL KEYS converge according to the convergence criteria defined in Main.py
        if self.n_iterations <= self.__start_iteration:
            has_converged = False
            print("Too soon to check convergence. I will iterate again.")
        else:
            abort = False  # Changes to True when next lines find a zone for a key that has not converged yet
                           # As long as there is ONLY ONE zone, code has not converged, so stop the loops as soon as you find that zone
            if self.__criteria == 'Previous':
                for zone in range(0,self.__n_zones):
                    for key in range(0,len(self.__names)):
                        difference = abs(self.__convResults[key][self.n_iterations][zone][0] - self.__convResults[key][self.n_iterations-1][zone][0])
                        average = 0.5 * (self.__convResults[key][self.n_iterations][zone][0] + self.__convResults[key][self.n_iterations-1][zone][0])
                        if not (difference < self.__tolerances[key]*average):
                            abort = True
                            break #break key loop
                    if abort: break #break zone loop

            elif self.__criteria == 'Variance':
                for zone in range(0,self.__n_zones):
                    for key in range(0,len(self.__names)):
                        mean = self.__convResults[key][self.n_iterations][zone][1]
                        variance = self.__convResults[key][self.n_iterations][zone][2]
                        delta_mean = np.sqrt(variance/(self.n_iterations-1))
                        if not (delta_mean < self.__tolerances[key]*mean):
                            abort = True
                            break
                    if abort: break
            elif self.__criteria == 'Median':
                for zone in range(0,self.__n_zones):
                    for key in range(0,len(self.__names)):
                        median = self.__convResults[key][self.n_iterations][zone][3]
                        delta_median = 0.5*self.__convResults[key][self.n_iterations][zone][4] /np.sqrt(self.n_iterations-1)
                        if not (delta_median < self.__tolerances[key]*median):
                            abort = True
                            break
                    if abort: break
            else:
                print("Warning: Convergence criterion not recognized! I will iterate until "+str(self.__max_iterations)+" have been done.")
            #Veredict
            if abort:
                has_converged = False
                print("Convergence not reached. I will iterate again.")
            else:
                has_converged = True
                print("Convergence achieved! Stopping...")

        #Increase number of iterations
        self.n_iterations += 1
        #And check if we have reached max iterations
        if self.n_iterations > self.__max_iterations:
            print("Warning: Stopping because max iterations reached.")
            print("I am doing this to avoid an infinite loop.")
            print("Check if tolerance is too low for the montecarlo noise.")
            print("If you want to improve resolution, try to:")
            print("-Increase the number of photon packets.")
            print("-Bias the probability distribution in those wavelength ranges you have more noise.")
            print("-Do more iterations to converge due to central limit theorem.")
            print("-Check statistics running post_process/averages.py .")
            has_converged = True

        return has_converged

    # TESTING

    def __diagnosticsInTable(self,filename):
        output = open(filename, 'w')
        for z in range(0,self.__n_zones):
            output.write("##### ZONE "+str(z)+" #####\n")

            output.write("iteration: ")
            for ii in range(self.n_iterations,0,-1):
                output.write(str(ii-1)+" ")
            output.write("\n")

            for key in range(0, len(self.__names)):
                output.write("## " + self.__names[key] + " ## \n")
                for type_of_data in range(0,5):
                    if type_of_data == 0:
                        output.write("convergence data: ")
                    elif type_of_data == 1:
                        output.write("means: ")
                    elif type_of_data == 2:
                        output.write("variances: ")
                    elif type_of_data == 3:
                        output.write("medians: ")
                    elif type_of_data == 4:
                        output.write("x84-x16 percentile error: ")
                    for ii in range(self.n_iterations-1, -1, -1):
                            output.write(str(self.__convResults[key, ii, z, type_of_data]) + " ")

                    output.write("\n")

        output.write("##########################################\n")
        output.write("That's all!")
        output.close()

    def showDiagnostics(self,filename='convergenceDiagnostics.dat'):
        self.__diagnosticsInTable(filename) #If you want a table

