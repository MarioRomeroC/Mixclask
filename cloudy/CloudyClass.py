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
from multiprocessing import Process
import numpy as np

class CloudyObject(converter.CloudyToSkirt):
    def __init__(self,options_dict):
        self.__initEsential(options_dict['AccuracyAndSpeed']['n_cpus'])
        self.__initChemistry()
        self.__initWavelengths(options_dict['Wavelength'])
        self.__initDetails()
        self.__parseData(options_dict['FileParameters']['ISM'])
        self.__ExePath = options_dict['Technical']['cloudy_path']
        self.__defineConstants() #Needed for ConvertedMethods
        
        self.__checkIssues()
        #self.showParams() #used for debug
    
    ### INITIALIZATION
    
    def __initEsential(self,n_cpus = 1):
        
        #Mandatory
        self._sedFiles = []
        self._param_mass = []
        self._param_nH = []
        self._n_zones = None
        
        #Geometry-dependend (not all of them are filled by run)
        self._geometry = [] #Each element of this list contains [Geometry-Type,params]
        #^^For example, [['ring',3,0.5,0.2],['shell',3.5,4.5]] indicates that first zone is a ring with their params, and second zone a shell with their params

        #Number of cpus
        self._n_cpus = n_cpus if n_cpus>1 else 1
        
    def __initChemistry(self):
        '''
        There are three layers of complexity:
        (1)    If helium abundance is not set (self._param_Y[zone] = None),
               I know that you are using the predefined cloudy options, where 
               you don't need any abundance.
        (2)    If metals abundance is set, He must be given as well, so
               I know that you aren't giving specific metals abundances.
        (3)    If metals abundances are set to None, but He is given, you are
               using specific abundances.
        As you noted, giving the metallicity have preference over custom abundances.
        '''
        
        #'metallicity' needs He (mass) fraction and metallicity to work
        self._param_X   = [] #H fraction (is computed from other params here)
        self._param_Y   = [] #He fraction
        self._param_Z   = [] #Metalicity
        
        #'custom' opens all elements modelled in cloudy (from He to Zn)
        #'Symbol':{'name':name,'A':mass number,'defaultToH':N_i/N_H,'abundance':mass_fraction}
        #'defaultToH' is the element fraction compared to H from table 7.4 of Hazy 1, it is used for default values
        self._param_element = {
            'H':{'name':'hydrogen','A':1.01,'defaultToH':1.0,'abundance': self._param_X},
            'He':{'name':'helium','A':4.00,'defaultToH':8.51e-2,'abundance':self._param_Y},
            'Li':{'name':'lithium','A':6.96,'defaultToH':1.12e-11,'abundance':[]},
            'Be':{'name':'beryllium','A':9.01,'defaultToH':2.40e-11,'abundance':[]},
            'B':{'name':'boron','A':10.8,'defaultToH':5.01e-10,'abundance':[]},
            'C':{'name':'carbon','A':12.0,'defaultToH':2.69e-4,'abundance':[]},
            'N':{'name':'nitrogen','A':14.0,'defaultToH':6.77e-5,'abundance':[]},
            'O':{'name':'oxygen','A':16.0,'defaultToH':4.90e-4,'abundance':[]},
            'F':{'name':'flourine','A':19.0,'defaultToH':3.63e-8,'abundance':[]},
            'Ne':{'name':'neon','A':20.2,'defaultToH':8.51e-5,'abundance':[]},
            'Na':{'name':'sodium','A':23.0,'defaultToH':1.74e-6,'abundance':[]},
            'Mg':{'name':'magnesium','A':24.3,'defaultToH':3.98e-5,'abundance':[]},
            'Al':{'name':'aluminium','A':27.0,'defaultToH':2.82e-6,'abundance':[]},
            'Si':{'name':'silicon','A':28.1,'defaultToH':3.24e-5,'abundance':[]},
            'P':{'name':'phosphorus','A':31.0,'defaultToH':2.57e-7,'abundance':[]},
            'S':{'name':'sulphur','A':32.1,'defaultToH':1.32e-5,'abundance':[]}, #'surfur' is the accepted name, but cloudy uses 'ph'
            'Cl':{'name':'chlorine','A':35.5,'defaultToH':3.16e-7,'abundance':[]},
            'Ar':{'name':'argon','A':40.0,'defaultToH':2.51e-6,'abundance':[]},
            'K':{'name':'potassium','A':39.1,'defaultToH':1.07e-7,'abundance':[]},
            'Ca':{'name':'calcium','A':40.1,'defaultToH':2.19e-6,'abundance':[]},
            'Sc':{'name':'scandium','A':45.0,'defaultToH':1.41e-9,'abundance':[]},
            'Ti':{'name':'titanium','A':47.9,'defaultToH':8.91e-8,'abundance':[]},
            'V':{'name':'vanadium','A':50.9,'defaultToH':8.51e-9,'abundance':[]},
            'Cr':{'name':'chromium','A':52.0,'defaultToH':4.37e-7,'abundance':[]},
            'Mn':{'name':'manganese','A':54.3,'defaultToH':2.69e-7,'abundance':[]},
            'Fe':{'name':'iron','A':55.8,'defaultToH':3.16e-5,'abundance':[]},
            'Co':{'name':'cobalt','A':58.9,'defaultToH':9.77e-8,'abundance':[]},
            'Ni':{'name':'nickel','A':58.7,'defaultToH':1.66e-6,'abundance':[]},
            'Cu':{'name':'copper','A':63.5,'defaultToH':1.55e-8,'abundance':[]},
            'Zn':{'name':'zinc','A':65.4,'defaultToH':3.63e-8,'abundance':[]},
            'Z':{'name':'metals','abundance':self._param_Z}
            }
        
        
        #Dust parameters
        self._param_DTG  = [] #Dust-to-gas
        self._param_qpah = [] #Pah-to-dust 
        
   
    def __initWavelengths(self,wavelength_dict):
        self._wavelength_max  = wavelength_dict['maxWavelength']
        self._wavelength_min  = wavelength_dict['minWavelength']
        self._wavelength_norm = wavelength_dict['normalization']
        self._wavelength_res  = wavelength_dict['resolution']
    
    def __massNumbersArray(self):
        result = []
        for symbol in self._param_element:
            if symbol != 'Z':
                A = self._param_element[symbol]['A']
                result.append(A)
        return np.array(result)
    
    def __defineConstants(self):
        #Mass number of elements from H to Zn
        self.Ai   = self.__massNumbersArray()
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
        #Default parameters when metallicity is given
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
        #Default options when specific (metal) abundances are given
        self._cloudy_fill_notGivenMetals = True #If true, elements not given will use 'abundances gass' value, and overwrite the values given
        
        #PAH options
        self._enable_PAH = True #Enables pah for the pre-defined mode
        self._disable_qheat = False #If true, you won't see the usual pattern that pah make in the spectrum, but cloudy is somewhat more stable.
        self._cloudy_default_qpah = 0.003962761126248864 #This is M_PAH/M_dust, or Pah-to-gas/dust-to-gas, assuming the same cloudy run as above.
        self._forced_qpah = 0.07 #Cloudy default value of qpah is too low compared with observational measures, which give about an order of magnitude higher (see figure 9 of Galliano+17) 
                                #Although this will be expanded as future work.
        #Cosmic rays
        self._use_cosmic_rays_background = True #If false and cloudy encounters molecular gas, it may crash
        #CMB
        self._use_CMB = True #Added for stability reasons. It will not propagate to Skirt, so you will not see the CMB in the output!
        
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
        
        #Extra outputs (They are used for diagnostics)
        self._ExtraCloudyOutputs = {   
            'Emissivity': False, #Stores emissivity j at last zone of cloudy
            'Opacity':False, #Stores opacity alpha at last zone of cloudy
            'RadiativeTransfer':{
                'Active': False, #Stores j, alpha and albedo as function of depth
                                        #for particular wavelengths defined below
                'Wavelengths': [200.0,443.0,550.0,2200,24e3,150e3] #nm
                },
           'Grains':{
               'Abundances': False, #Save grains abundances (g/cm3) of all species (dust and pah) as function of depth
               'DTG':False #Save the dust-to-gas ratio of all species as function of depth
       }
        
    }

    ### MAKE AND RUN CLOUDY INPUT

    def make_input(self,outfile="cloudyInput"):
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
            self.__writeExtra(outfile,z,self._ExtraCloudyOutputs)
            
            outfile.close()
            self.__inputs.append(filename.replace('.in','')) #Will be needed for running cloudy
    
    def run_input(self):
        #Old code
        '''
        for zone in range(0,len(self.__inputs)):
            currInput = self.__inputs[zone]
            os.system(self.__ExePath+" -r "+currInput)
            if self.__problemDisaster(currInput):
                self.__errorHandleInput(zone,currInput)
            #else do nothing
        '''
        def execute_input(input,region):
            os.system(self.__ExePath + " -r " + input)
            if self.__problemDisaster(input):
                self.__errorHandleInput(region, input)
            # else do nothing
        n_inputs = len(self.__inputs)
        n_cpus = self._n_cpus
        for zone in range(0,n_inputs,n_cpus):
            currPrograms = []
            for s in range(0,n_cpus):
                try: currInput = self.__inputs[zone+s]
                except IndexError: break #No more inputs
                p = Process(target=execute_input,args=(currInput,zone+s))
                currPrograms.append(p)
                p.start()

            #Do not start next iteration until above processes finish!
            for program in currPrograms:
                program.join()

    
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
                if self._param_qpah[zone] != None and self._param_qpah[zone] > 0.0:
                    grains_scale *= (1.0-self._param_qpah[zone])
                file.write(str(grains_scale))
                if no_qheat:
                    file.write(" no qheat \n")
                else:
                    file.write(" \n")
                
                #PAH
                if self._param_qpah[zone] != None and self._param_qpah[zone] > 0.0: #Repeated for legibility
                    file.write("grains pah ")
                    pah_scale = self._param_qpah[zone]/self._cloudy_default_qpah
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
            #You only gave M and n_H. No info for other chemical elements
            file.write("abundaces gass \n")
            
        elif self._param_Z[zone] != None:
            #metallicity is given AND TAKES PREFERENCE over custom abundances
            #first, use solar abundances (see table 7.4 of cloudy hazy 1)
            #THESE ABUNDANCES DOES NOT ADD GRAINS (see table 7.2 of cloudy hazy 1)
            file.write("abundances gass \n")
            #He content
            scale_He = self._param_Y[zone] / self._cloudy_default_Y
            file.write("element helium scale "+str(scale_He)+" \n")
            #Metal content
            scale_Z = self._param_Z[zone]/self._cloudy_default_Z
            file.write("metals "+str(scale_Z)+" \n")
        
        else: #You gave custom abundances
            if self._cloudy_fill_notGivenMetals:
                file.write("abundances gass \n") #Elements given will overwrite these values
            else:
                file.write("abundances all -20 \n") #Elements NOT given will be negligible
            #Now set each element
            H_mass_fraction = self._param_element['H']['abundance'][zone]
            for symbol in self._param_element:
                elem_name = self._param_element[symbol]['name']
                elem_mass_fraction = self._param_element[symbol]['abundance'][zone]
                if symbol!='H' and elem_mass_fraction != None and symbol!='Z':
                    mass_number = self._param_element[symbol]['A']
                    #None = not given
                    elem_relative_to_H = elem_mass_fraction / ( H_mass_fraction*mass_number )
                    value = np.log10(elem_relative_to_H)
                    file.write("element "+elem_name+" abundance "+str(value)+" \n")
                
        
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
        
        if self._use_CMB:
            outfile.write("CMB \n")
        
        sedfile.close()
    
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
    #Chemistry parameters
    def __fillElemen(self,symbol,data):
        self._param_element[symbol]['abundance'].append(float(data) if data != None else None)
    
    #Dust parameters (optional)
    def __fillDTG(self,data):
        self._param_DTG.append(float(data) if data != None else None)
    def __fillqPah(self,data):
        self._param_qpah.append(float(data) if data != None else None)
    
    
    def __fillUnfilled(self):
        #Chemical elements
        for symbol in self._param_element:
            if len(self._param_element[symbol]['abundance']) == 0:
                for zone in range(0,self._n_zones):
                    self.__fillElemen(symbol, None)
        #Dust
        if len(self._param_DTG) == 0:
            for zone in range(0,self._n_zones):
                self.__fillDTG(None)
                self.__fillqPah(None)
    
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
            #Chemistry-Basic (if some parameters are lacking, it will use default options)
            'helium': self.__fillHe,
            'metallicity': self.__fillZ,
            #Chemistry-Custom
            'lithium':'Li',
            'beryllium':'Be',
            'boron':'B',
            'carbon':'C',
            'nitrogen':'N',
            'oxygen':'O',
            'flourine':'F',
            'neon':'Ne',
            'sodium':'Na',
            'magnesium':'Mg',
            'aluminium':'Al',
            'silicon':'Si',
            'phosphorus':'P',
            'sulphur':'S',
            'sulfur':'S', #Both words are correct for this element
            'chorine':'Cl',
            'argon':'Ar',
            'potassium':'K',
            'calcium':'Ca',
            'scandium':'Sc',
            'titanium':'Ti',
            'vanadium':'V',
            'chromium':'Cr',
            'manganese':'Mn',
            'iron':'Fe',
            'cobalt':'Co',
            'nickel':'Ni',
            'copper':'Cu',
            'zinc':'Zn',
            #Dust (optional, it not given, it'll take 'grains ism' option)
            #Give '0' or negative number to disable dust
            'dusttogas': self.__fillDTG,
            'dust-to-gas': self.__fillDTG,
            'pah-to-dust': self.__fillqPah,
            'pahtodust': self.__fillqPah,
            'q_pah': self.__fillqPah
            }
        #Do
        if isinstance(switch[key],str):
            self.__fillElemen(switch[key], data)
        else: #I'm a function
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
        
        #And compute Hydrogen fraction
        for i in range(0,self._n_zones):
            self._param_element['H']['abundance'][i] = self._computeHydrogenFraction(i)
        
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
    def __checkIssues(self):
        #This method revises if something odd is found.
        for zone in range(0,self._n_zones):
            if self._param_DTG[zone] != None and self._param_DTG[zone] < 0.0:
                print("Warning: Negative dust-to-gas ratio in zone "+str(zone)+". Is this intended? I will ignore dust there")
            if self._param_qpah[zone] != None and self._param_qpah[zone] < 0.0:
                print("Warning: Negative pah-to-dust ratio (q_pah) in zone "+str(zone)+". Is this intended? I will ignore dust there")
            for symbol in self._param_element:
                if symbol == 'H' or symbol == 'He' or symbol == 'Z': 
                    if self._param_element[symbol]['abundance'][zone] != None and self._param_element[symbol]['abundance'][zone] < 0.0:
                        raise RuntimeError("Error: Negative "+symbol+" abundance. I cannot go on, revise your input.")
                else:
                    if self._param_element[symbol]['abundance'][zone] != None and self._param_element[symbol]['abundance'][zone] < 0.0:
                        print("Warning: Negative "+symbol+" abundance in zone "+str(zone)+". Is this intended? I will ignore this abundance there")
                        #I haven't write a elem_abundance <0.0 in the pertinent method, so I correct here
                        self._param_element[symbol]['abundance'][zone] = None
    
    
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
            
        
        
    
