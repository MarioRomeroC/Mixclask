#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Auxiliary functions to move and copy documents

Mario Romero            July 2021
'''

import os
import numpy as np

# =============================================================================
# ROUTINES WITH TERMINAL
# =============================================================================

move = lambda file,des : os.system("mv "+file+" -t "+des)
copy = lambda file,des : os.system("cp "+file+" -t "+des)
runPython = lambda script,argv='' : os.system("python3 "+script+argv)
changedir = lambda des: os.system("cd "+des)

def generateSkirtFile(scriptFolder,scriptFile):
    changedir(scriptFolder)
    runPython(scriptFile)
    os.system("cd -") #Return to original directory

# =============================================================================
# ROUTINES WITH LISTS
# =============================================================================

def relist(my_str):
    #Geometry data comes in lists : ['type', param1, param2, etc]
    #When read from a file, the list is my_str = "['type', param1, param2, etc]"
    #This function is to revert this
    
    #Remove some characters
    tmp_list = my_str.replace('[','')
    tmp_list = tmp_list.replace(']','')
    tmp_list = tmp_list.replace("'",'')
    
    #Create the list
    tmp_list = tmp_list.split(',')
    
    #Convert in floats any number inside my_list
    my_list = []
    for elem in tmp_list:
        try:
            my_list.append(float(elem))
        except ValueError:
            my_list.append(elem)
    
    return my_list

# =============================================================================
# ROUTINES WITH ARRAYS
# =============================================================================

mean = lambda array : 0.5*(array[1:] + array[:-1])
def downsample(x,X,Y,interp_type='Integral'): 
    #Note that len(y) = len(x)-1. use x=mean(x) later if you want x to have same length
    ### Generate interpolation ###
    
    X = np.array(X)
    Y = np.array(Y)
    y = []
    
    if interp_type != 'Integral':
        #Do the interpolation
        y = np.interp(x,X,Y)
    else:
        #Compute the integral F = int(Y*dX)
        dX = np.diff(X)
        X  = mean(X)
        Y  = mean(Y)
        F  = np.cumsum(Y*dX)
        #Compute the interpolation OF THE INTEGRAL, f
        f = np.interp(x,X,F)
        #Compute the derivative of f to get y
        y = np.diff(f)/np.diff(x)
        #x = mean(x)
        #print(len(f),len(x),len(y))
    return y


def testDownsample(x,y,X,Y):
    #make sure that len(x)=len(y)
    import matplotlib.pyplot as plt
    plt.plot(X,Y,'r-')
    plt.plot(x,y,'b-')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

takeFirst = lambda String : String.split()[0]
#def takeFirst(String):
#    string_data = String.split()
#    return string_data[0]

def RMSD(old,new):
    #Root-Mean-Squared-Deviation of two arrays
    #Make sure that len(old) = len(new)
    x = np.array(old)
    y = np.array(new)
    
    RMSD = np.sqrt((sum((y-x)*(y-x)))/len(x) )
    return RMSD

def integrate(x_lim,x,y):
    #Integrate according to trapezoidal rule, len(x) = len(y)
    #Unlike np.trapz, x is NOT the whole region of integration
    x0 = x_lim[0]
    y0 = np.interp(x0,x,y)
    xf = x_lim[1]
    yf = np.interp(xf,x,y)
    
    suma = 0.0
    for i in range(1,len(x)):
        if x[i-1] <= x0 and x0 < x[i]:
            #Initial step
            suma += 0.5*(x[i]-x0)*(y[i]+y0)
        elif x0 <= x[i-1] and x[i] < xf:
            suma += 0.5*(x[i]-x[i-1])*(y[i]+y[i-1])
        elif x[i-1] <= xf and xf < x[i]:
            #Last step
            suma += 0.5*(xf-x[i-1])*(yf-y[i-1])
    
    return suma

# =============================================================================
# ROUTINES WITH FILES
# =============================================================================    

def readColumn(filename,column):
    '''
    This returns the column selected (in python notation, so column 1 is '0' 
    and so on) from the file.
    column can be a list if you want to extract multiple columns, although
    data will have a different format.
    '''
    file = open(filename,'r')
    
    data = []
    while True:
        line = file.readline()
        if not line: break #EoF
        elif line[0] == '#': continue #commented line
        
        line_data = line.split()
        if not isinstance(column,list):
            data.append(float(line_data[column]))
        else:
            row = [float(line_data[col]) for col in column]
            data.append(row)
    
    file.close()
    return np.array(data)

def addColumn(matrix,column):
    #Here you add column to matrix.
    
    result = []
    for i in range(0,len(column)):
        row = np.append(matrix[i] , column[i])
        result.append(row)
    return np.array(result)

def findNormalization(x_value,x,y):
    '''
    This return an array of [x_norm, y(x_norm)]
    x_norm is the value of x[i] that x[i+1] > x_value but x[i]<=x_value
    y can be a matrix, where each column is a different function
    '''
    #return np.append(x_value,np.interp(x_value,x,y))
    for i in range(0,len(x)-1):
        #print(x[i+1],x_value,y[:,i+1])
        if (x[i+1] >= x_value) and (x[i] < x_value):
            try:
                return np.append(x[i+1],y[:,i+1])
            except IndexError:
                #This route is taken if y is 1d (you give a list, not a 'matrix' as commented above)
                return np.array([x[i+1],y[i+1]])
    raise RuntimeError("Normalization not found!")

def createProbabilityFile(original_file_path,original_folder,probability_folder='',invert=False):
    #This creates a probability file from 'original_file'.
    #For skirt reasons, they need a file with two columns: One for wavelengths, other for its probability density.
    #Probability density does not need to be normalized, but it needs to be given in luminosity units or skirt crashes.

    original_file = original_file_path.replace(original_folder+'/','')
    output_name = original_file.replace('.stab','_probability.stab')
    output_name = output_name.replace('.txt', '_probability.stab')
    output_name = output_name.replace('.dat', '_probability.stab')

    # Create file header
    output = open(output_name,'w')
    output.write("# Column 1: wavelength (nm) \n")
    output.write("# Column 2: probability density (erg/s) \n")
    output.write("# luminosity units are required for skirt to work. Probability density has no units \n")
    # Extract columns
    columns = readColumn(original_file_path, [0, 1])
    wavelength = columns[:,0]
    density = columns[:,1]
    if invert:
        density = 1.0/density

    #Write file
    for l in range(0,len(wavelength)):
        output.write(str(wavelength[l])+" "+str(density[l])+" \n")

    output.close()
    #Move output to probability file
    move(output_name,probability_folder)
    return probability_folder+'/'+output_name
