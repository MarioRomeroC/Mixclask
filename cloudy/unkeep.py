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

#DEPRECIATED!
def GenerateCloudyFiles(Jfilename):
    '''
    This script reads the file generated with 'MeanIntensityAtPositions' probe,
    and generates a readable input for cloudy
    '''
    from cloudy_params import CloudyParams
    
    n_digits = 4
    '''
    Known bug: If you give a table with too many significant digits to cloudy,
    it may crash due to overshoot in cloudy interpolation routine 
    (cloudy rounds/truncates all numbers from the table, and it may found that
     some results are higher than the interpolation limits and crash).
    You can avoid round your results by giving a higher number to 'n_digits'
    and hope that cloudy does not crash randomly in one particular zone (as it
    may happen).
    '''
    
    useIntensity = True #If True, it will use the cloudy intensity at range
                            #Otherwise, it will use nuf(nu) X at Y (Ryd)
    wl_nuF = 1e7 #nm (used if above is false)
    
    microns_TO_nm = 1000.
    W_per_m2_TO_erg_per_cm2 = 1000.
    R_inf = 10973731.57
    nm_to_Ryd = lambda x : 1.0/(x*1e-9*R_inf)
    #Read file
    file = open(Jfilename,'r')
    wavelength = []
    x = []
    y = []
    z = []
    nuJnu = []
    while True:
        line = file.readline()
        if not line: #EoF
            break
        line_data = line.split()
        
        if line_data[0] == '#':
            #Here we colect the wavelengths.
            #The structure of the skirt output is
            #Column 1: lambda*J_lambda at lambda = 0.0001 micron (W/m2/sr)
            #If we split line in line_data, the only element that can be
            #converted to a float is 0.0001, which is the wavelength
            for element in line_data:
                try:
                    wavelength.append(float(element))
                except:
                    continue
        else:
            nuJnu.append([])
            #First 3 columns are x,y,z
            x.append(float(line_data[0]))
            y.append(float(line_data[1]))
            z.append(float(line_data[2]))
            for i in range(3,len(line_data)):
                nuJnu[-1].append(float(line_data[i]))
    file.close()
    
    x = np.array(x)
    y = np.array(y)
    z = np.array(z)
    #Also flip next arrays, as cloudy wants it in opposite order
    #wavelength = np.flip(np.array(wavelength)) * microns_TO_nm #in nm
    #nuJnu = np.fliplr(np.array(nuJnu))
    wavelength = np.array(wavelength) * microns_TO_nm #in nm
    nuJnu = np.array(nuJnu) * W_per_m2_TO_erg_per_cm2
    
    #Cloudy may fail if you give too much significant digits
    round_to = lambda n,x : round(x, -int(np.floor(np.log10(abs(x)))) + (n - 1))
    
    #Create the sed files
    outputfiles = CloudyParams['files'] 
    #make sure that the order is respected, file 1 is for position 1 and so on
    '''
    The number of outputs are tuned in the number of lines in 'MeanIntensity_Positions.txt'.
    If you do not want problems, they MUST match both the number of zones and the order:
        file 1 is for first UNCOMMENTED (i.e: no #) line 1, and so on
    Extra outputs for J may be useful as well, and they are the following lines AFTER THE ZONES.
    They will be labelled as 'extra_Xpc.sed','extra_Ypc.sed' and so on
    '''
    if len(nuJnu) < len(outputfiles):
        raise RuntimeError("Number of outputs does not match!")
    for i in range(0,len(nuJnu)):
        R = np.sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i])
        fluxfile = None
        if i < len(outputfiles):
            fluxfile = open(outputfiles[i],'w')
        else:
            fluxfile = open("extra_"+str(R)+"pc.sed",'w')
        #fluxfile.write("# nuFnu at "+str(R)+" pc \n")
        fluxfile.write("# column 1: wavelength (nm) \n")
        fluxfile.write("# column 2: 4pi*nu*J_nu (erg/cm2/s) \n")
        if useIntensity:
            #Control_value is the integral at wavelengths of intensity
            total_J = np.trapz(nuJnu[i]/wavelength,wavelength)
            fluxfile.write("# intensity "+str(round_to(n_digits,np.log10(4.0*np.pi*total_J)))+" range "+str(wavelength[0])+" to "+str(wavelength[-1])+" nm \n")
        else:
            #Control value is nu*F at wl_nuF
            fixed_values = findNormalization(wl_nuF,wavelength,nuJnu)
            #Above line may make coudy to fail. If that happens, cut the decimals in the respective cloudy input and repeat
            #fixed_values = np.append(wl_nuF,np.interp(wl_nuF,wavelength,nuJnu))
            #print(fixed_values)
            fixed_nuF = 4.0*np.pi*fixed_values[i+1] #Check if that 4pi is needed...
            photon_energy = nm_to_Ryd(fixed_values[0])
            fluxfile.write("# nuf(nu) "+str(round_to(n_digits,np.log10(fixed_nuF)))+" at "+str(photon_energy)+" \n")
            #Cloudy hazy says that nuF is actually nuJ.
            
        options_written = False
        for j in range(len(nuJnu[i])-1,0,-1):
            nuJnu_write = nuJnu[i][j]
            if nuJnu_write < 1e-300: nuJnu_write = 1e-300 #Cloudy cannot read 0.0
            fluxfile.write(str(round_to(n_digits,wavelength[j]))+" "+str(round_to(n_digits,4.0*np.pi*nuJnu_write)))
            if not options_written:
                fluxfile.write(" nuFnu units nm \n")
                options_written = True
            else:
                fluxfile.write(" \n")
        fluxfile.close()
    