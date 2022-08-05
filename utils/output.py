#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This script allows to generate a radiation field that is result of several iterations, computing its average.

Mario Romero        August 2022
'''

import unkeep as unk
from glob import glob
import os
import numpy as np
np.seterr(all='raise')

###### FUNCTION MODE ######

def compute_average(pointers,folder='',i_0=1,i_f=None):

    wavelengths = None
    all_individual_data = []
    header = ''
    for ii in range(i_0,i_f):
        it_folder= folder+'iteration'+str(ii)+'/'
        it_individual_data = []
        for file in pointers:
            file_path = it_folder + file
            #Collect header
            if header == '':
                infile = open(file_path,'r')
                rows = 0 #I only expect to fill 2 rows
                while True:
                    line = infile.readline()
                    if 'column' in line.lower() and rows<2:
                        header += line
                        rows += 1
                    else: break
                infile.close()
            #Collect data
            data = unk.readColumn(file_path,[0,1])
            wavelengths  = data[:,0]
            it_individual_data.append( data[:,1] )
        all_individual_data.append(it_individual_data)
    all_individual_data = np.array(all_individual_data)
    #To make indexing easier for you and I
    it_to_index = lambda iteration: iteration - i_0

    #Create output files
    for z in range(0,len(pointers)):
        output = open('Average_'+pointers[z],'w')
        output.write(header)
        output.write("# column 3: relative error () \n")
        for l in range(0,len(wavelengths)):
            output.write(str(wavelengths[l])+" ")
            #compute mean
            mean = sum(all_individual_data[:,z,l]) / len(all_individual_data[:,z,l])
            output.write(str(mean)+" ")

            #compute relative error
            variance = sum( (all_individual_data[:,z,l] - mean )**2.0 )/ (i_f-i_0-1.0)
            error = np.sqrt(variance) / mean
            output.write(str(error)+" ")

            output.write("\n")
        output.close()
    #move output to folder
    os.system("mv Average_* -t "+folder)



folder_path = '../' #This is the path to the folder that contains all iterations (it has several subfolders called 'iteration0','iteration1',...
#You need to know the filenames you want to make an average
#You can give a list, or use glob (e.g.: glob(folder_path+'*.sed') will do the trick in most cases)
filenames   = glob(folder_path+'*.sed')
#If you have NOT used glob as suggested, comment next line
filenames = [filenames[i].replace(folder_path,'') for i in range(0,len(filenames))]
start_iteration = 2 #Do not consider iterations below this number
end_iteration   = 5 #Do not consider iterations above this number
#Also, be aware that all 'filenames' in all iterations must have the same wavelength grid (this is not an issue if it comes from the same Mixclask run)

compute_average(filenames,folder_path,start_iteration,end_iteration)