#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This script allows to generate a radiation field that is result of several iterations.
It computes its average and other statistical data.

Mario Romero        August 2022
'''

#from ..utils.unkeep import readColumn
from glob import glob
import os
import numpy as np
np.seterr(all='warn')

#Set your locations
inputs_path = '../YourResultsFolder/'
outputs_folder = 'Statistics' #It will be created in inputs_path
outputs_key = 'Stat'
#Set starting and ending iteration
i_0 = 1
i_f = np.infty #If absurd numbers, this code will use all iterations that are found in inputs_path
#Set percentiles (if desired)
percentiles = [16,84] #Median is always given, no need to add a '50' in this list


########## OTHER FUNCTIONS ##########
def extract_number(string):
    result = ''
    for char in string:
        if char.isdigit():
            result += char
    return int(result)

########## ROUTINE ##########
#find folders and files
folders = sorted(glob(inputs_path+"iteration*"),key=extract_number)
#^They are sorted by iteration (iteration0 is element 0, iteration1 element 1, and so on)
filenames = glob(folders[0]+'/'+'*.sed')
for f in range(0,len(filenames)):
    filenames[f] = filenames[f].replace(folders[0]+'/','')
    #Collect data
    wavelength = None #You only need to extract this data once
    wl_extracted = False
    all_iterations_data = []
    for ii in range(0,len(folders)):
        #Skip not desired data
        iteration = extract_number(folders[ii])
        if iteration<i_0 or iteration>i_f: continue

        #Get proper data
        if not wl_extracted:
            wavelength, FourPi_nuJ = np.loadtxt(folders[ii]+'/'+filenames[f],unpack=True,skiprows=2,usecols=(0,1))
            wl_extracted = True
        else:
            FourPi_nuJ = np.loadtxt(folders[ii]+'/'+filenames[f],unpack=True,skiprows=2,usecols=(1))

        all_iterations_data.append(FourPi_nuJ)
    all_iterations_data = np.array(all_iterations_data)

    #Create outputs
    output = open(outputs_key+'_'+filenames[f],'w')
    output.write("# column 1: wavelength (nm) \n")
    output.write("# statistics of 4pi*nu*J_nu (erg/cm2/s) \n")
    output.write("# column 2: mean \n")
    output.write("# column 3: variance \n")
    output.write("# column 4: median \n")
    col = 5
    try:
        for per in percentiles:
            output.write("# column "+str(col)+": percentile "+str(per)+" \n")
            col += 1
    except IndexError:
        pass #Do nothing
    output.write("# Number of iterations : "+str(len(all_iterations_data))+" \n")

    for wl in range(0,len(wavelength)):
        #Create statistical data
        mean   = np.mean(all_iterations_data[:,wl])
        var    = np.var(all_iterations_data[:,wl])
        median = np.median(all_iterations_data[:,wl])
        output.write(str(wavelength[wl])+" "+str(mean)+" "+str(var)+" "+str(median)+" ")
        try:
            for per in percentiles:
                percentile = np.percentile(all_iterations_data[:,wl],per)
                output.write(str(percentile)+" ")
        except IndexError:
            pass #Do nothing
        output.write("\n")

    output.close()

#Create results folder
os.system("mkdir "+inputs_path+outputs_folder)
#And move all there
os.system("mv "+outputs_key+"* -t "+inputs_path+outputs_folder)
#That's all!
