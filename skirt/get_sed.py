#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 17:17:45 2021

@author: pablo
"""
from glob import glob
from numpy import sort


def get_from_folder(path,files=None):
    #Mario: I'm adding a default variable for backwards compatibility 
    #(following functions use THIS function)
    #If files = None, it will use the original behaviour
    #Else files = [file1,file2,etc] so paths = [path/file1,path/file2,etc]
    if files is None:
        paths = sort(glob(path+'/*')) 
    else:
        paths = [path+'/'+files[i] for i in range(0,len(files))]
    return paths

def get_norm(path, units, files=None):
    paths = get_from_folder(path,files)
    all_norms = []
    for path_i in paths:
        with open(path_i, 'r') as f:
                all_lines = f.readlines()                                
                normalization = all_lines[3].split(' ')[2]
                normalization += ' '+units
        f.close()
        all_norms.append(normalization)
    return all_norms
def get_norm_wavelength(path, units, files=None):
    paths = get_from_folder(path,files)
    all_wave_norms = []
    for path_i in paths:
        with open(path_i, 'r') as f:
                all_lines = f.readlines()                
                norm_wavelength = all_lines[2].split(' ')[3]
                norm_wavelength += ' '+units
        f.close()
        all_wave_norms.append(norm_wavelength)
    return all_wave_norms

def get_optdepth_wavelength(path, units, files=None):
    paths = get_from_folder(path,files)
    all_wave_norms = []
    for path_i in paths:
        with open(path_i, 'r') as f:
                all_lines = f.readlines()                
                norm_wavelength = all_lines[4].split(' ')[3]
                norm_wavelength += ' '+units
        f.close()
        all_wave_norms.append(norm_wavelength)
    return all_wave_norms
                
def get_optdepth_norm(path, files=None):
    paths = get_from_folder(path,files)
    all_norms = []
    for path_i in paths:
        with open(path_i, 'r') as f:
                all_lines = f.readlines()                
                norm = float(all_lines[5].split(' ')[2])
        f.close()
        all_norms.append(norm)
    return all_norms                       
                
#Mario: I add this function to allow mass normalization
def get_mass_norm(path, units, files=None):
    paths = get_from_folder(path,files)
    all_mass_norms = []
    for path_i in paths:
        with open(path_i, 'r') as f:
                all_lines = f.readlines()                
                norm_wavelength = all_lines[4].split(' ')[3]
                norm_wavelength += ' '+units
        f.close()
        all_mass_norms.append(norm_wavelength)
    return all_mass_norms

