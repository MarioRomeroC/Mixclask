#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This class generates a cloudy input file, and also executes them.
I'm taking lines from Pablo's script

Mario Romero            July 2021
'''

import os
#import glob
from utils.unkeep import relist
import cloudy.ConverterMethods as converter
#import src.ParserClass as parser
import numpy as np

class CloudyObject(converter.CloudyToSkirt):
    def __init__(self,params_path,cloudy_path,wavelength_dict):
        self.__initParams()
        self.__initWavelengths(wavelength_dict)
        self.__parseData(params_path)
        self.__ExePath = cloudy_path
        self.__defineConstants() #Needed for ConvertedMethods
        self.__initDetails()
        
        #self.showParams() #used for debug
    
    ### INITIALIZATION
    
    def __initParams(self):
        
        #Mandatory
        self._sedFiles = []
        self._param_mass = []
        self._param_nH = []
        self._n_zones = None
        
        #Chemistry and dust
        self._param_Y   = [] #He fraction
        self._param_Z   = [] #Metalicity
        self._param_DTG = [] #Dust-to-gas
        
        #Geometry-dependend (not all of them are filled by run)
        self._geometry = [] #Each element of this list contains [Geometry-Type,params]
        #^^For example, [['ring',3,0.5,0.2],['shell',3.5,4.5]] indicates that first zone is a ring with their params, and second zone a shell with their params
        
        
        
    def __initWavelengths(self,wavelength_dict):
        self._wavelength_max  = wavelength_dict['maxWavelength']
        self._wavelength_min  = wavelength_dict['minWavelength']
        self._wavelength_norm = wavelength_dict['normalization']
        self._wavelength_res  = wavelength_dict['resolution']
    
    def __defineConstants(self):
        #Mass number of elements from H to Zn
        self.Ai   = np.array([1.0, 4.0, 6.9, 9.0, 11., 12., 14., 16., 19., 20., 23., 24., 27., 28., 31., 32., 35., 40., 39., 40., 45., 48., 51., 52., 55., 56., 59., 59., 64., 65.])
        #Physical constants
        self.mH   = 1.67353e-24    #g [Hydrogen mass]
        self.pc   = 3.08568e18     #cm
        self.Msun = 1.988e33       #g [Solar mass]
        self.Rinf = 10973731.57   #m-1 [Rydberg constant]
        #Relevant conversions
        self.microns_TO_nm = 1000.
        self.W_per_m2_TO_erg_per_cm2 = 1000.
        self.nm_to_Ryd = lambda x : 1.0/(x*1e-9*self.Rinf)
        
    def __initDetails(self):
        #This option is used for debug purposes only
        self._use_intensity_command_in_cloudy = False #This will use intensity X at range instead of nuf(nu) Y at Z
        #Default parameters
        '''
        These parameters are taken with a cloudy run with these commands:
            abundances gass
            grains ism
            grains pah
            set pah constant
            table ism
            stop zone 1
        Other options shown in this method are taken as FALSE and thus not included
        '''
        self._cloudy_default_Y   = 0.25057390087996967 #He abundance by default
        self._cloudy_default_Z   = 0.0134294404311891 #metals abundance by default
        self._cloudy_default_DTG = 0.006606 #dust-to-gas ratio scales with this number. Includes grains and pah
        
        #PAH options
        self._enable_PAH = True
        self._disable_qheat = False #If true, you won't see the usual pattern that pah make in the spectrum, but cloudy is somewhat more stable.
        self._cloudy_default_qpah = 0.003962761126248864 #This is M_PAH/M_dust, or Pah-to-gas/dust-to-gas, assuming the same cloudy run as above.
        self._forced_qpah = 0.05 #Default value of qpah is too low compared with observational measures, which give about an order of magnitude higher (see figure 9 of Galliano+17) 
                                #Although this will be expanded as future work.
        #Cosmic rays
        self._use_cosmic_rays_background = True #If false and cloudy encounters molecular gas, it may crash
        #CMB
        self._use_CMB = False
        
        '''
        Known bug: If you give a table with too many significant digits to cloudy,
        it may crash due to overshoot in cloudy interpolation routine 
        (cloudy rounds/truncates all numbers from the table, and it may found that
         some results are higher than the interpolation limits and crash).
        You can avoid round your results by giving a higher number to 'n_digits'
        and hope that cloudy does not crash randomly in one particular zone (as it
        may happen).
        '''
        self._n_digits = 4
        self.round_to = lambda n,x : round(x, -int(np.floor(np.log10(abs(x)))) + (n - 1))

    ### MAKE AND RUN CLOUDY INPUT

    def make_input(self,IOextensions,outfile="cloudyInput"):
        #Generates the input file for cloudy
        #IOextensions is a dictionary
        prefix = outfile
        self.__inputs = []
        for z in range(0,self._n_zones):
            sed_filename  = self._sedFiles[z]
            
            filename = prefix+str(z)+".in"
            
            outfile = open(filename,'w')
            outfile.write("title cloudy zone "+str(z)+" \n")
            self.__writeChemistry(outfile,z,self._disable_qheat)
            self.__writeRadiation(outfile, sed_filename)
            self.__writeGeometry(outfile,z)
            self.__writeOptions(outfile)
            self.__writeOutputs(outfile,z) #Needed outputs from cloudy
            
            #Extra outputs?
            self.__writeExtra(outfile,z,IOextensions)
            
            outfile.close()
            self.__inputs.append(filename.replace('.in','')) #Will be needed for running cloudy
    
    def run_input(self):
        for zone in range(0,len(self.__inputs)):
            currInput = self.__inputs[zone]
            os.system(self.__ExePath+" -r "+currInput)
            if self.__problemDisaster(currInput):
                self.__errorHandleInput(zone,currInput)
            #else do nothing
    
    def __writeChemistry(self,file,zone,no_qheat):
        #Mandatory parameters
        hden = np.log10(self._param_nH[zone]) #cloudy prefers the logarithm
        
        def __writeDust():
            #Encapsulated here, as here is the only place where you use it
            if self._param_DTG[zone] == None: #Default option
                file.write("grains ism ")
                if no_qheat:
                    file.write("no qheat \n")
                else:
                    file.write("\n")
                #Add pah, if desired
                if self._enable_PAH:
                    file.write("grains pah ")
                    if self._disable_qheat:
                        file.write("no qheat \n")
                    else:
                        file.write("\n")
                    file.write("set pah constant \n")
            elif self._param_DTG[zone] > 0.0:
                '''
                When dust is included in this model, Z is assumed to be for gas only.
                If pah is enabled, dust is splitted in grains and pah according to q_PAH (see 'init_details()' for the numbers).
                The cloudy input should look like this:
                    grains ism DTG*(1-qpah)
                    grains pah DTG*qpah
                DTG is (DTG_given/self._default_DTG). (cloudy always asks for a number to rescale from its default, and does not want an absolute number)
                In 'grains ism', qpah = self._forced_qpah.
                In 'grains pah', qpah = (self._forced_qpah/self._default_qpah). As pah need to be rescaled from the default value of cloudy (see 'init_details')
                Therefore, the DTG you give will be the DTG reported in cloudy with 'save grain D/G ratio'
                '''
                file.write("grains ism ") #It's the default, 'grains' does the same
                grains_scale = self._param_DTG[zone]/self._cloudy_default_DTG #This is the DTG from the '''
                if self._enable_PAH:
                    grains_scale *= (1.0-self._forced_qpah)
                file.write(str(grains_scale))
                if no_qheat:
                    file.write(" no qheat \n")
                else:
                    file.write(" \n")
                
                #PAH
                if self._enable_PAH: #Repeated for legibility
                    file.write("grains pah ")
                    pah_scale = self._forced_qpah/self._cloudy_default_qpah
                    pah_scale *= self._param_DTG[zone]/self._cloudy_default_DTG
                    file.write(str(pah_scale))
                    if no_qheat:
                        print("Warning: the characteristic PAH emission won't be reflected if qheat is off")
                        file.write(" no qheat \n")
                    else:
                        file.write(" \n")
                    #To make this value constant across a cloudy region
                    file.write("set pah constant \n")
                
            #else do nothing
                
        file.write("## CHEMICAL COMPOSITION \n")
        file.write("hden "+str(hden)+" \n")
        
        if self._param_Y[zone] == None: #Default options
            file.write("abundaces ism no grains \n")
            __writeDust()
            
        else:
            #if Z is given but not specific metals
            #first, use solar abundances (see table 7.4 of cloudy hazy 1)
            #THESE ABUNDANCES DOES NOT ADD GRAINS (see table 7.2 of cloudy hazy 1)
            file.write("abundances gass \n")
            #He content
            scale_He = self._param_Y[zone] / self._cloudy_default_Y
            file.write("element helium scale "+str(scale_He)+" \n")
            #Metal content
            scale_Z = self._param_Z[zone]/self._cloudy_default_Z
            file.write("metals "+str(scale_Z)+" \n")
            
            __writeDust()
            
        
        
        if self._use_cosmic_rays_background:
            file.write("cosmic rays background \n")
            
    def __writeRadiation(self,outfile,sedname):
        outfile.write("## RADIATION FIELD \n")
        outfile.write('''table sed "'''+sedname+'''" \n''')
        #Open sedfile and get the intensity
        sedfile = open(sedname,'r')
        while True:
            line = sedfile.readline()
            if line[0] != '#' or not line:
                break
            elif 'intensity' in line:
                line = line.replace('# ', '')
                outfile.write(line)
                break
            elif 'nuf(nu)' in line:
                line = line.replace('# ', '')
                outfile.write(line)
                break
        sedfile.close()
        
        if self._use_CMB:
            outfile.write("CMB \n")
    
    def __writeGeometry(self,file,zone):
        file.write("## GEOMETRY \n")
        thickness = self._getThickness(zone) #In ConvertedMethods
        
        file.write("stop thickness "+str(thickness)+" parsec linear \n")
    
    def __writeOptions(self,file):
        file.write("## OTHER OPTIONS \n")
        #file.write("iterate to convergence \n") #USED FOR DEBUGGING
        file.write("stop temperature off \n")
    
    def __writeOutputs(self,file,zone,output_ext=".txt"):
        file.write("## OUTPUTS \n")
        file.write('''save overview "overview_zone'''+str(zone)+str(output_ext)+'''" last \n''')
        file.write('''save continuum "spectra_zone'''+str(zone)+str(output_ext)+'''" last units nm \n''')
        file.write('''save abundances "composition_zone'''+str(zone)+str(output_ext)+'''" last \n''')
        file.write('''save optical depth "tau_zone'''+str(zone)+str(output_ext)+'''" last units nm \n''')
    
    def __writeExtra(self,file,zone,my_dict):
        output_ext = ".txt"
        if my_dict['Emissivity']:
            file.write('''save diffuse continuum "emissivity_zone'''+str(zone)+str(output_ext)+'''" last units nm \n''')
        if my_dict['Opacity']:
            file.write('''save total opacities "opacity_zone'''+str(zone)+str(output_ext)+'''" last units nm \n''')
        if my_dict['RadiativeTransfer']['Active']:
            for wavelength in my_dict['Wavelengths']:
                file.write('''save continuum emissivity '''+str(self.nm_to_Ryd(wavelength)))
                file.write(''' "radTransfer_'''+str(wavelength)+"nm_zone"+str(zone)+str(output_ext)+'''" last \n''')
        if my_dict['Grains']['Abundances']:
            file.write('''save grain abundance "grainAbundance_zone'''+str(zone)+str(output_ext)+'''" last units nm \n''')
        if my_dict['Grains']['DTG']:
            file.write('''save grain D/G ratio "grainDTG_zone'''+str(zone)+str(output_ext)+'''" last units nm \n''')
    
    ### PARSE DATA ### 
    def __fillGeometry(self,data):
        #data = ['type',params]
        #Because data would be a string and I want a list, something else needs to be done first
        true_data = relist(data)
        self._geometry.append(true_data)
    def __fillSED(self,data):
        #data = filename
        self._sedFiles.append(data)
    def __fillMass(self,data):
        #data = float
        self._param_mass.append(float(data))
    def __fillnH(self,data):
        self._param_nH.append(float(data))
    #Optional parameters
    def __fillHe(self,data):
        self._param_Y.append(float(data) if data != None else None)
    def __fillZ(self,data):
        self._param_Z.append(float(data) if data != None else None)
    def __fillDTG(self,data):
        self._param_DTG.append(float(data) if data != None else None)
    
    
    def __fillUnfilled(self):
        if len(self._param_Y) == 0:
            for zone in range(0,self._n_zones):
                self.__fillHe(None)
        if len(self._param_Z) == 0:
            for zone in range(0,self._n_zones):
                self.__fillZ(None)
        if len(self._param_DTG) == 0:
            for zone in range(0,self._n_zones):
                self.__fillDTG(None)
    
    def __fillData(self,key,data):
        #This expands the list self._header_order
        
        switch = {
            #Geometry
            'geometry': self.__fillGeometry,
            #Output
            'sedfile': self.__fillSED,
            'outputfile': self.__fillSED,
            #Chemistry-Mandatory
            'mass': self.__fillMass,
            'hydrogendensity': self.__fillnH,
            #Chemistry-Other (if some parameters are lacking, it will use default options)
            'hellium': self.__fillHe,
            'metallicity': self.__fillZ,
            #Dust (optional, it not given, it'll take 'grains ism' option)
            #Give '0' or negative number to disable dust
            'dusttogas': self.__fillDTG,
            'dust-to-gas': self.__fillDTG
            }
        #Do
        switch[key](data)
    
    def __parseData(self,file_path):
        file = open(file_path,'r')
        header = []
        
        #Read the file
        while True:
            line = file.readline()
            if not line: break #EoF
            line_data = line.split()
            
            if line_data[0] == '#': 
                #I should expect: '# Column X : Key (units)
                header.append(line_data[4].lower())
                
            else:
                #I should expect: data1 data2 data3 data4 ...
                #Get relevant data
                for i in range(0,len(header)):
                    self.__fillData(header[i], line_data[i])
                    
              
        #And the final param
        self._n_zones = len(self._sedFiles)
        #Unfilled parameters remain as '[]', giving len = 0. 
        #Nevertheless, I will fill with 'None', to make handling with these options easier
        self.__fillUnfilled()
        
        file.close()
    
    ### GIVE SOMETHING ###
    
    def giveSEDfiles(self):
        return self._sedFiles
    
    ### ERROR HANDLING ###       
    
    def __problemDisaster(self,filename):
        # Checks if 'PROBLEM DISASTER' (phrase that appears when cloudy crashes) appears in the output
        outfile = open(filename+'.out','r')
        
        disaster_found = False
        while True:
            line = outfile.readline()
            if not line:
                break #EoF
            elif 'PROBLEM DISASTER' in line:
                disaster_found = True
                break #cloudy crashed, no need for more reading
        
        outfile.close()
        return disaster_found
    
    def __errorHandleInput(self,zone,filename):
        #This method is used to avoid this script for crashing when cloudy crashes (if possible).
        #Something, the crash is avoid by removing some options by the input.
        
        if self._disable_qheat == False:
            '''
            qheat error:
                self._disable_qheat = False enables quantum heating in cloudy.
                Sometimes, the cloudy code crashes when running related processes.
                This is solved by disabling qheat with 'no qheat' in the cloudy inputs
            '''
            #Rewrite the input
            os.system("rm "+filename+".in")
            outfile = open(filename+".in",'w')
            outfile.write("title cloudy zone "+str(zone)+" no qheat run \n")
            self.__writeChemistry(outfile,zone,True) #I'm disabling qheat here
            self.__writeRadiation(outfile, self._sedFiles[zone])
            self.__writeGeometry(outfile,zone)
            self.__writeOptions(outfile)
            self.__writeOutputs(outfile,zone) #Needed outputs from cloudy
            #No extra outputs
            outfile.write("# Extra outputs are disabled for second, error handling, runs")
            outfile.close()
            
            #Rerun input
            os.system(self.__ExePath+" -r "+filename)
            #And recheck if problem is solved
            if self.__problemDisaster(filename):
                #Not solved
                raise RuntimeError("Cloudy crashed in zone "+str(zone)+". Removing qheat did not solve the issue. Check "+filename+".out for more details.")
            else:
                #Solved
                print("Warning: Cloudy crashed in zone "+str(zone)+". I removed qheat in this zone and solved the issue.")
            
        else:
            raise RuntimeError("Cloudy crashed in zone "+str(zone)+". Check "+filename+".out for more details.")
    
    ### DEBUG ###
    
    def showParams(self):
        print("-Number of zones : "+str(self._n_zones))
        
        print("Wavelength parameters:")
        print("-Wavelength range : "+str(self._wavelength_min)+" to "+str(self._wavelength_max)+" nm")
        print("-Wavelength resolution : "+str(self._wavelength_res))
        print("-Wavelength of reference : "+str(self._wavelength_norm)+" nm \n")
        
        print("Each zone data:")
        for i in range(0,self._n_zones):
            print("---ZONE "+str(i)+"---")
            geometry_type = self._geometry[i][0]
            print("Geometry : "+geometry_type)
            if geometry_type == 'shell':
                print("-Inner radius : "+str(self._geometry[i][1])+" pc")
                print("-Outer radius : "+str(self._geometry[i][2])+" pc")
            elif geometry_type == 'ring':
                print("-Ring radius : "+str(self._geometry[i][1])+" pc")
                print("-Ring width  : "+str(self._geometry[i][2])+" pc")
                print("-Ring height : "+str(self._geometry[i][3])+" pc")
            
            print("Hydrogen density : "+str(self._param_nH[i])+" cm-3")
            print("Total mass       : "+str(self._param_mass[i])+" Msun")
            Y = self._param_Y[i]
            if Y != None:
                print("Hellium fraction : "+str(Y))
            else:
                print("Hellium fraction : Not specified")
            Z = self._param_Z[i]
            if Y != None:
                print("Metallicity      : "+str(Z))
            else:
                print("Metallicity      : Not specified")
            
            DTG = self._param_DTG[i]
            if DTG != None and DTG > 0.0:
                print("Dust-to-gas      : "+str(DTG))
            elif DTG == None:
                print("Dust-to-gas      : Not specified")
            else:
                print("Dust-to-gas      : Not included")
            
        
        
    
