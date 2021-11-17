#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import cloudy.CloudyClass as cc #class with all methods needed to operate with cloudy
from skirt import pigs #source code that generates the .ski file
from skirt.ski_params import SkiParams #class that contains all the data to generate the .ski file
import os
import time
from cloudy.unkeep import move #utility functions

# Where are the gas and star parameters starting from this folder?
gas_params  = 'input_data/params/Your_gas_file.dat'
star_params = 'input_data/params/Your_star_file.dat'  
meanIntensity_positions = 'input_data/MeanIntensity_Positions/Your_positions.txt'
# Do we start from the beginning?
is_iteration0 = True # If true, mixclask will do a skirt simulation with only 'star_params' data before running cloudy
                     # If false, you should have the results of the radiation field of a previous run in root folder (same as this file)
                     #      For the latter case, it is useful if you want to do more iterations than originally intended.

# What are your spectra resolution and normalization?
WavelengthOptions = {
    'normalization': 1.0e7, #nm
    'maxWavelength':1.0e9, #nm
    'minWavelength':0.1, #nm
    'resolution':200
    }

#These outputs are useful for diagnostics, they are not needed for this code to work properly
ExtraCloudyOutputs = {   
    'Emissivity': True, #Stores emissivity j at last zone of cloudy
    'Opacity':True, #Stores opacity alpha at last zone of cloudy
    'RadiativeTransfer':{
        'Active': False, #Stores j, alpha and albedo as function of depth
                                #for particular wavelengths defined below
        'Wavelengths': [200.0,443.0,550.0,2200,24e3,150e3] #nm
        },
   'Grains':{
       'Abundances': True, #Save grains abundances (g/cm3) of all species (dust and pah) as function of depth
       'DTG':True #Save the dust-to-gas ratio of all species as function of depth
       }
        
    }

# Some technical parameters
cloudy_path = '/path/to/your/cloudy/exe'
show_cloudy_params = False
last_iteration = 0 #Useful for tracking the number of iteration if 'is_iteration0=False'
n_iterations = 3
n_threads = 4

### MAIN ROUTINE ###

t_start = time.time()
cloudy = cc.CloudyObject(gas_params,cloudy_path,WavelengthOptions)
skirt_params = SkiParams(gas_params,star_params,meanIntensity_positions,WavelengthOptions)
#cloudy.showParams()
#raise
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
for it in range(1,n_iterations+1):
    # =============================================================================
    # Create cloudy
    # =============================================================================
    
    folder = "iteration"+str(it+last_iteration)
    
    t_it = time.time()
    print("iteration "+str(it+last_iteration)+" ("+str(t_it-t_start)+" s):")
    
    
    cloudy.make_input(ExtraCloudyOutputs)
    
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
    
    
    
    
# iteration0 -> while[skirt -> conversion1 -> cloudy -> conversion2 -> check convergence -> skirt]
if show_cloudy_params: cloudy.showParams()
t_end = time.time()
print("Program finished! ("+str(t_end-t_start)+" s) \n")
