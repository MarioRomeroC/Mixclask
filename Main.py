#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import cloudy.CloudyClass as cc #class with all methods needed to operate with cloudy
from skirt import pigs #source code that generates the .ski file
from skirt.ski_params import SkiParams #class that contains all the data to generate the .ski file
import utils.convergence as conv
from utils.unkeep import move #utility functions
import os
import time

# Where are the gas and star parameters starting from this folder?
gas_params  = 'input_data/params/Your_gas_file.dat'
star_params = 'input_data/params/Your_star_file.dat'  
meanIntensity_positions = 'input_data/MeanIntensity_Positions/Your_positions.txt'
# Do we start from the beginning?
is_iteration0 = True # If true, mixclask will do a skirt simulation with only 'star_params' data before running cloudy
                     # If false, you should have the results of the radiation field of a previous run in root folder (same as this file)
                     #  For the latter case, it is useful if you want to do more iterations than originally intended.

# What are your spectra resolution and normalization?
WavelengthOptions = {
    'normalization': 550.0, #nm
    'maxWavelength':3.0e5, #nm
    'minWavelength':10.0, #nm
    'resolution':200, #Equispaced in log
    'convWavelength':[(10.0,90.0),(100e3,300e3)] #nm #Convergence wavelength
    }

# Some technical parameters
cloudy_path = '/path/to/your/cloudy/exe'
show_cloudy_params = False
n_iterations = 15
n_threads = 2
tolerance = 0.15 #For convergence.

### MAIN ROUTINE ###

t_start = time.time()
cloudy = cc.CloudyObject(gas_params,cloudy_path,WavelengthOptions)
skirt_params = SkiParams(gas_params,star_params,meanIntensity_positions,WavelengthOptions)
program = conv.ConvergenceObject(cloudy.giveSEDfiles(),WavelengthOptions['convWavelength'],n_iterations,tolerance)

#cloudy.showParams()

### ITERATION 0 ###
if is_iteration0:
    last_iteration = 0
    t_it = time.time()
    print("iteration 0 ("+str(t_it-t_start)+" s):")
    skirt_params.prepareSkiFile(is_iteration0)
    pigs.SkirtFile(skirt_params, output_path='skirt_file')
    
    os.system("skirt -t "+str(n_threads)+" skirt_file.ski > tmp.txt") #Make sure you followed skirt instructions
    cloudy.GenerateCloudyFiles("skirt_file_nuJnu_J.dat")
    folder = "iteration0"
    os.system("mkdir "+folder)
    os.system("mv *.dat -t "+folder)
    os.system("cp *.sed -t "+folder) #copy because I need these files in the root folder
    os.system("mv *.ski -t "+folder)

### FOLLOWING ITERATIONS ###
t_it = time.time()
#for it in range(1,n_iterations+1):
while(not program.stop()):
    # =============================================================================
    # Create cloudy
    # =============================================================================
    
    folder = "iteration"+str(program.iteration(t_it-t_start))
    
    cloudy.make_input()
    
    cloudy.run_input()
    
    cloudy.GenerateSkirtInput()
    
    #Save the data in a separate folder    
    print("Moving cloudy data to "+folder)
    os.system("mkdir "+folder)
    os.system("mv *.txt -t "+folder)
    os.system("cp *.stab -t "+folder) #copy and then move
    os.system("mv *.in -t "+folder)
    os.system("mv *.out -t"+folder)
    
    #Move the relevant data to 'input_data' folder
    move("MeanFileGasMix_*.stab","input_data/gas_props")
    move("GasSource_*.stab","input_data/gas_sources")
    
    # =============================================================================
    # Create skirt
    # =============================================================================
    
    
    skirt_params.prepareSkiFile()
    pigs.SkirtFile(skirt_params, output_path='skirt_file')
    
    #print("skirt -t "+str(n_threads)+" skirt_file.ski > tmp.txt")
    os.system("skirt -t "+str(n_threads)+" skirt_file.ski > tmp.txt") #Make sure you followed skirt instructions
    
    cloudy.GenerateCloudyFiles("skirt_file_nuJnu_J.dat")
    
    
    print("Moving skirt data to "+folder+" \n")
    os.system("mv *.dat -t "+folder)
    os.system("cp *.sed -t "+folder) #copy because I need these files in the root folder
    os.system("mv *.ski -t "+folder)
    
    t_it = time.time()
    
    
# iteration0 -> while[skirt -> conversion1 -> cloudy -> conversion2 -> check convergence -> skirt]
if show_cloudy_params: cloudy.showParams()
t_end = time.time()
print("Program finished! ("+str(t_end-t_start)+" s) \n")
