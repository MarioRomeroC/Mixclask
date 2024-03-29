#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 25 11:12:18 2021

@author: pablo
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from skirt.get_sed import get_from_folder, get_norm, get_norm_wavelength, get_optdepth_norm, get_optdepth_wavelength, get_mass_norm
from flatten_dict import flatten, unflatten
from utils.unkeep import relist,createProbabilityFile
import numpy as np

def iterate_over_dict(dictionary):
    for key, value in dictionary.items():
        if isinstance(value, dict):
            for pair in iterate_over_dict(value):
                yield (key, *pair)
        else:
            # If value is not dict type then yield the value
            yield (key, value)

def translate_dictionary(my_dict,ii):
        new_dict = dict()
        for prop_jj in iterate_over_dict(my_dict):
            keywords = tuple( [prop_jj[elem] for elem in range(0,len(prop_jj)-1)] )
            value = prop_jj[-1]
            if isinstance(value,(list,np.ndarray)):
                value = value[ii]
            elif value is None:
                value = {}
            new_dict[keywords] = value
        new_dict = unflatten(new_dict)
        dict_type = list(new_dict.keys())[0]
        return [new_dict[dict_type],dict_type]

              
# =============================================================================
# PARAMETER FILE FOR MAKING SKIRT .SKI FILES
# =============================================================================
class SkiParams(object):
    def __init__(self,options_dict):#gas_path,star_path,output_positions,wavelength_dict):
        self.__initParams()
        self.__initWavelengths(options_dict['Wavelength'])
        self.__initFolders(options_dict['FileParameters']['stars']['folder'])
        self.__initPhotons(options_dict['AccuracyAndSpeed']['photon_packets'],options_dict['AccuracyAndSpeed']['PhotonProbability'])
        self.__initDetails()
        
        self.__parseData(options_dict['FileParameters']['ISM'],'Gas')
        self.__parseData(options_dict['FileParameters']['stars']['file'],'Star')
        
        self.__deduceLimits()
        self.__deduceLinks()
        self.__checkSkirtIssues()
        self.__Positions = options_dict['FileParameters']['positions']
    
    ### INITIALIZATION
    
    def __initParams(self):
        #Geometry GAS params
        self._gas_geometry  = []
        #'shell' geometry = ['shell',Rin,Rout] -> Rin = inner radius, Rout = outer radius
        #'ring' geometry = ['ring',R,w,h] -> R=middle radius, w=width, h=height
        #number of gas zones
        self._gas_zones = None

        #Geometry STAR params
        self._star_files = []
        self._star_geometry = []
        #'shell' geometry -> Same as gas
        #'ring' geometry -> Same as gas
        #'point' geometry -> ['point',x,y,z]
        #number of stellar zones
        self._star_zones = None
        # Are stellar zones related to gas zones?
        self._stars_gas_link = []
        # > Index of this array is the stellar zone
        # > Number inside each value is the index of the gas zone, 'None' is used if stellar region does not have gas counterpart.
        # > This is to check if some stars are inside gas, for resolution purposes.
    
    def __initFolders(self,folder_path):
        #folder locations
        # DO NOT CHANGE THESE PARAMETERS, THEY ARE HARDCODED BECAUSE CLOUDY MOVES THE FILES THERE
        self._gas_opacity_folder  = "input_data/gas_props"
        self._gas_sources_folder  = "input_data/gas_sources"
        # Except next one. As you see, you can define its path in Main.py
        self._star_sources_folder = folder_path #"input_data/star_sources"
    
    def __initWavelengths(self,wavelength_dict):
        self._wavelength_max  = wavelength_dict['maxWavelength']
        self._wavelength_min  = wavelength_dict['minWavelength']
        self._wavelength_norm = wavelength_dict['normalization']
        self._wavelength_res  = wavelength_dict['resolution']

    def __initPhotons(self,n_photons,probability_dic):
        self._photonPackets = n_photons
        self._perRegionProbability = probability_dic['per_region']
        self._distributionBias = probability_dic['wavelengthBias']
        self._customProbabilityFile = probability_dic['customDistributionFile']
        #Details
        self._probability_folder = 'input_data/probability_distributions'

    def __initDetails(self):
        #This function contains parameters used for debugging
        self._montecarloSeed = '0'
        self._fluxProbeDistance = '1 Mpc'
        self._geometryUnits = ' pc' #The space IS intended
        self._resolution_r = 200
        self._resolution_z = 200
        #Skirt cries if you use 'ring' geometry and place one at R = 0pc
        self._skirt_ringRadius_correction = 0.01 #pc
        self._z_ring_limit = 10.0 #In case of all gas-zones are 'ring' geometry, The upper limit would be this number multiplied by the maximum height (h) found in all zones 
        #Iteration0 options
        self._nullMaterialFile = 'input_data/params/NullMaterial.stab' #Skirt will use this file as 'material' in the iteration0 (should be a material without absorption and scattering)
        self._nullMass = 1.0e-10 #Msun. Normalization of 'NullMaterial'. The ideal value is 0.0, but I selected a small number because skirt will not run instead
    
    def __deduceLimits(self):
        #limits are decided with gas zones.
        #We first need to find the greater radius
        R_max = 0.0
        z_max = 0.0
        
        for ii in range(0,self._gas_zones):
            Rcan = None
            if self._gas_geometry[ii][0] == 'shell':
                Rcan = self._gas_geometry[ii][2] #Rout
                zcan = Rcan #Spherical symmetry!
            elif self._gas_geometry[ii][0] == 'ring':
                Rcan = self._gas_geometry[ii][1] + self._gas_geometry[ii][2] #R+w
                zcan = self._gas_geometry[ii][3]*self._z_ring_limit #h*self._z_ring_limit
            else:
                raise RuntimeError("Geometry not implemented!")
            
            if Rcan > R_max:
                R_max = Rcan
            if zcan > z_max:
                z_max = zcan
        
        for ii in range(0,self._star_zones):
            Rcan = None
            if self._star_geometry[ii][0] == 'shell':
                Rcan = self._star_geometry[ii][2] #Rout
                zcan = Rcan #Spherical symmetry!
            elif self._star_geometry[ii][0] == 'ring':
                Rcan = self._star_geometry[ii][1] + self._gas_geometry[ii][2] #R+w
                zcan = self._star_geometry[ii][3]*self._z_ring_limit #h*self._z_ring_limit
            elif self._star_geometry[ii][0] == 'point':
                x = self._star_geometry[ii][1]
                y = self._star_geometry[ii][2]
                z = self._star_geometry[ii][3]

                if any((x,y,z)) != 0.0:
                    raise RuntimeError("Sorry, I have not implemented points outside the origin =( ")
                    #It would require a change in the coordinate system used (Cylindrical coordinates are no longer valid)
                    #And to allow that is now low priority

                Rcan = np.sqrt(x*x+y*y+z*z)
                zcan = z
        
            if Rcan > R_max:
                R_max = Rcan
            if zcan > z_max:
                z_max = zcan
        
        self._border_r = str(R_max)+self._geometryUnits
        self._border_z = str(z_max)+self._geometryUnits

    def __deduceLinks(self):
        #This function checks if stellar zones are related with gas zones.
        for ii in range(0,self._star_zones): #Per stellar zone
            if self._star_geometry[ii][0] == 'point':
                x = self._star_geometry[ii][1]
                y = self._star_geometry[ii][2]
                z = self._star_geometry[ii][3]
                #Check if it is inside a gas zone
                found_link = False
                for jj in range(0,self._gas_zones): #Per gas zone
                    #Compute parameters
                    r, Rmin, Rmax = None, None, None
                    if self._gas_geometry[jj][0] == 'shell':
                        r    = np.sqrt(x*x+y*y+z*z)
                        Rmin = self._gas_geometry[jj][1]
                        Rmax = self._gas_geometry[jj][2]
                    elif self._gas_geometry[jj][0] == 'ring':
                        r    = np.sqrt(x*x+y*y)
                        Rmin = self._gas_geometry[jj][1]-self._gas_geometry[jj][2]
                        Rmax = self._gas_geometry[jj][1]+self._gas_geometry[jj][2]
                        h    = self._gas_geometry[jj][3]
                        #false if point is above ring. Skip final if
                        if abs(z) > h : continue
                    else: raise RuntimeError("Geometry not implemented")
                    #Check if related
                    if r <= Rmax and r >= Rmin:
                        self._stars_gas_link.append(jj)
                        found_link = True
                        break
            elif self._star_geometry[ii][0] == 'shell':
                for jj in range(0,self._gas_zones):
                    Rmin, Rmax, W = None, None, None
                    r = (self._star_geometry[ii][1] + self._star_geometry[ii][2]) / 2.0
                    w = (self._star_geometry[ii][2] - self._star_geometry[ii][1])
                    if self._gas_geometry[jj][0] == 'shell':
                        Rmin = self._gas_geometry[jj][1]
                        Rmax = self._gas_geometry[jj][2]
                        W    = self._gas_geometry[jj][2] - self._gas_geometry[jj][1]
                    elif self._gas_geometry[jj][0] == 'ring':
                        continue #A ring is a very small region compared with a shell! There is no point!
                    else: raise RuntimeError("Geometry not implemented")
                    #Check if related
                    if r <= Rmax and r >= Rmin and w <= W:
                        self._stars_gas_link.append(jj)
                        found_link = True
                        break
            elif self._star_geometry[ii][0] == 'ring':
                Rmin, Rmax, W = None, None, None
                r = self._star_geometry[ii][1]
                w = self._star_geometry[ii][2]
                h = self._star_geometry[ii][3]
                for jj in range(0,self._gas_zones):
                    if self._gas_geometry[jj][0] == 'shell':
                        Rmin = self._gas_geometry[jj][1]
                        Rmax = self._gas_geometry[jj][2]
                        W = self._gas_geometry[jj][2] - self._gas_geometry[jj][1]
                    elif self._gas_geometry[jj][0] == 'ring':
                        W = self._gas_geometry[jj][2]
                        Rmin = self._gas_geometry[jj][1]-W
                        Rmax = self._gas_geometry[jj][1]+W
                    #Check if related
                    if r <= Rmax and r >= Rmin and w <= W:
                        self._stars_gas_link.append(jj)
                        found_link = True
                        break
            if not found_link:
                self._stars_gas_link.append(None)

            #And that's it


    def __checkSkirtIssues(self):
        #(1) If ring geometry is set and some of 'ringRadius' is 0.0
        #   Skirt crashes due they do not accept rings with R=0
        radius_correction = 0.01 #self._geometry_units
        for ii in range(0,self._gas_zones):
            if self._gas_geometry[ii][0] == 'ring' and self._gas_geometry[ii][1] == 0.0:
                self._gas_geometry[ii][1] = radius_correction
        
        #Same issue with stars
        for ii in range(0,self._star_zones):
            if self._star_geometry[ii][0] == 'ring' and self._star_geometry[ii][1] == 0.0:
                self._star_geometry[ii][1] = radius_correction
    
    ### PARSING DATA
    
    def __parseData(self,file_path,Type):
        #This goes a little different than 'parseData' in CloudyClass.py
        #Type has two options: 'Gas' or 'Star'
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
                #If this is the star file
                if Type == 'Star':
                    self._star_files.append(line_data[header.index('sedfile')])
                    self._star_geometry.append(relist(line_data[header.index('geometry')]))
                elif Type == 'Gas':
                    #only this input is needed for skirt
                    self._gas_geometry.append(relist(line_data[header.index('geometry')]))
                
        file.close()
        
        #Now tell me if data is 'gas' or 'star'
        if Type == 'Gas':
            self._gas_zones  = len(self._gas_geometry)
        elif Type == 'Star':
            self._star_zones = len(self._star_geometry)
        
        #Now yes, the end!
    
    ### CREATE DICTIONARIES FOR SKIRT
    #Original code made by Pablo Corcho-Caballero
    
    def prepareSkiFile(self,iteration0=False):    
            
            #Mario: For legibility, I split this method in subrutines
            self.__createBasics()
            self.__createSources(iteration0)
            self.__createMedia(iteration0)
            self.__createInstruments()
            self.__createProbes()
        
    def __createBasics(self):
        # =============================================================================
        # BASIC PROPERTIES 
        # =============================================================================
        Basics = {
                    'MonteCarloSimulation':{
                        'userLevel':"Expert",
                        'simulationMode':"ExtinctionOnly",
                        'iterateMediumState':'false',
                        'iterateSecondaryEmission':'false',
                        'numPackets':self._photonPackets,
                        'random':{
                            'type':'Random',
                            'Random':{'seed':self._montecarloSeed}
                                    },
                        'units':{
                            'type':'Units',
                            'ExtragalacticUnits':{
                                'wavelengthOutputStyle': 'Wavelength',
                                'fluxOutputStyle': 'Neutral'
                                    }
                                },
                        'cosmology':{
                            'type':'Cosmology',
                            'LocalUniverseCosmology':{}
                                    }
                    }
                }
        #Return
        self.Basics = Basics

    def __wavelengthBiasDistribution_options(self,zone,are_stars=True,iteration0=False):
        result_dic = {
            'type' : 'WavelengthDistribution',
        }

        #Two possible options to write
        def writeSkirtDefault():
            subresult = {
                    'minWavelength': str(self._wavelength_min) + ' nm',
                    'maxWavelength': str(self._wavelength_max) + ' nm'
                }
            return subresult
        def writeFileDefault(filename):
            subresult = {
                'filename': filename
            }
            return subresult


        if self._perRegionProbability == 'Custom':
            result_dic['FileWavelengthDistribution'] = writeFileDefault(self._customProbabilityFile)
        # If user has not provided its own file, then I first look for cases in which the default logWavelength distribution from skirt is used
        elif self._perRegionProbability == 'logWavelength' or iteration0:
            #This fires if
            #1- 'logWavelength' is selected in main (AccuracyAndSpeed->PhotonProbability->per_region : 'logWavelength')
            #2- iteration 0

            result_dic['LogWavelengthDistribution'] = writeSkirtDefault()
        elif self._perRegionProbability == 'Luminosity':
            filename = ''
            # Straightforward. Just use the luminosity of each zone
            if are_stars: filename = (get_from_folder(self._star_sources_folder, self._star_files))[zone]
            else: filename = (get_from_folder(self._gas_sources_folder))[zone]

            result_dic['FileWavelengthDistribution'] = writeFileDefault(filename)
        elif self._perRegionProbability == 'invLuminosity':
            filename = ''
            folder = ''
            # Straightforward. Just use the luminosity of each zone
            if are_stars:
                folder = self._star_sources_folder
                filename = (get_from_folder(folder, self._star_files))[zone]
            else:
                folder = self._gas_sources_folder
                filename = (get_from_folder(folder))[zone]
            #And add this extra step
            filename = createProbabilityFile(filename,folder,self._probability_folder,invert=True)
            result_dic['FileWavelengthDistribution'] = writeFileDefault(filename)

        elif are_stars and self._stars_gas_link[zone] is None:
            #An exception of following options. If region is a stellar one and DO NOT have a gas counterpart, use 'logWavelength' in them
            result_dic['LogWavelengthDistribution'] = writeSkirtDefault()
        elif self._perRegionProbability == 'Extinction':
            # Find the right media_file
            folder = self._gas_opacity_folder
            true_zone = self._stars_gas_link[zone] if are_stars else zone #self._gas_gas_link[zone] is not None because it is the previous elif
            media_file = (get_from_folder(folder))[true_zone]

            #filename = media_file #PLACEHOLDER!
            filename = createProbabilityFile(media_file,folder,self._probability_folder)
            result_dic['FileWavelengthDistribution'] = writeFileDefault(filename)
        else:
            print("Warning: Probability distribution have not been recognized, using 'logWavelength' option")
            result_dic['LogWavelengthDistribution'] = writeSkirtDefault()

        return result_dic

    def __createSources(self,iteration0):
        # =============================================================================
        # SOURCES 
        # =============================================================================
        Sources = dict()
        
        # Basic properties 
        Sources['SourceSystem'] = {
                                    'minWavelength': str(self._wavelength_min)+' nm',
                                    'maxWavelength': str(self._wavelength_max)+' nm',
                                    'wavelengths': str(self._wavelength_norm)+' nm',
                                    'sourceBias': '0.5'
                                    }
        
        Sources['SourceSystem']['sources'] = {'type':'Source'}

        # STELLAR SOURCES
        n_skirt_sources_count = 0 #star+gas sources.
        
        starSource_properties = None
        all_star_files        = get_from_folder(self._star_sources_folder, self._star_files)
        stars_norm_wavelength = get_norm_wavelength(self._star_sources_folder, 'nm', self._star_files), #It will check the third line of 'fileSED', check if units are correct in your file
        stars_norm_value      = get_norm(self._star_sources_folder, 'erg/s', self._star_files) #It will check the fourth line of 'fileSED',check if units are correct in your file
        
        #print(stars_norm_wavelength)
        for ii in range(0,self._star_zones):
            n_skirt_sources_count += 1
            #print(ii)
            star_file_ii  = all_star_files[ii]
            norm_wl_ii    = [stars_norm_wavelength[0][ii]]
            norm_value_ii = [stars_norm_value[ii]]
        
            if self._star_geometry[ii][0] == 'shell':
                Rin  = str(self._star_geometry[ii][1])+self._geometryUnits #e.g.: '3.0 pc'
                Rout = str(self._star_geometry[ii][2])+self._geometryUnits
                
                starSource_properties = {'GeometricSource':{ #This line indicates the type of source it is. Inside this dictionary you have ALL parameters needed for that source
                    		                'velocityMagnitude':'0 km/s', 
                    		                'sourceWeight':"1", 
                    		                'wavelengthBias': self._distributionBias,
                    		                
                    		                'geometry':{
                    		                    'type':'Geometry',
                    		                    'ShellGeometry':{
                        				            'minRadius': [Rin],  
                        				            'maxRadius': [Rout], 
                        				            'exponent': ['0']
                        		                        }
                    		                    },
                    		                'sed':{
                    		                    'type':'SED',
                    		                    'FileSED':{'filename':[star_file_ii]}
                    		                    },
                    		                'normalization':{
                    		                    'type':'LuminosityNormalization',
                    		                    'SpecificLuminosityNormalization':{
                    		                            'wavelength':norm_wl_ii,
                    		                            'unitStyle':'neutralmonluminosity',
                    		                            'specificLuminosity':norm_value_ii} 
                    		                    },
                                            'wavelengthBiasDistribution':self.__wavelengthBiasDistribution_options(ii,iteration0=iteration0)
                    		                }
                    		            }
            elif self._star_geometry[ii][0] == 'ring':
                R = str(self._star_geometry[ii][1])+self._geometryUnits
                w = str(self._star_geometry[ii][2])+self._geometryUnits
                h = str(self._star_geometry[ii][3])+self._geometryUnits
                
                starSource_properties = {'GeometricSource':{ #This line indicates the type of source it is. Inside this dictionary you have ALL parameters needed for that source
                    		                'velocityMagnitude':'0 km/s', 
                    		                'sourceWeight':"1", 
                    		                'wavelengthBias': self._distributionBias,
                    		                
                    		                'geometry':{
                    		                    'type':'Geometry',
                    		                    'RingGeometry':{
                        				            'ringRadius': [R], 
                        				            'width': [w], 
                        				            'height': [h]
                    		                        }
                    		                    },
                    		                'sed':{
                    		                    'type':'SED',
                    		                    'FileSED':{'filename':[star_file_ii]}
                    		                    },
                    		                'normalization':{
                    		                    'type':'LuminosityNormalization',
                    		                    'SpecificLuminosityNormalization':{
                    		                            'wavelength':norm_wl_ii, 
                    		                            'unitStyle':'neutralmonluminosity',
                    		                            'specificLuminosity':norm_value_ii} 
                    		                    },
                                            'wavelengthBiasDistribution':self.__wavelengthBiasDistribution_options(ii,iteration0=iteration0)
                    		                }
                    		            }
                
            elif self._star_geometry[ii][0] == 'point':
                x = str(self._star_geometry[ii][1])+self._geometryUnits
                y = str(self._star_geometry[ii][2])+self._geometryUnits
                z = str(self._star_geometry[ii][3])+self._geometryUnits
                
                starSource_properties = {'PointSource':{
                                            'positionX': [x],
                                            'positionY': [y],
                                            'positionZ': [z],
                                            'velocityX': ['0 km/s'],
                                            'velocityY': ['0 km/s'],
                                            'velocityZ': ['0 km/s'],
                                            'sourceWeight': '1',
                                            'wavelengthBias': self._distributionBias,
                                            'angularDistribution': {
                                                'type':'AngularDistribution',
                                                'IsotropicAngularDistribution':None #This is to copy '{}' in the dictionary
                                                },
                                            'polarizationProfile': {
                                                'type': 'PolarizationProfile',
                                                'NoPolarizationProfile':None
                                                },
                                                
                                            'sed':{ 'type':'SED',
                                                'FileSED':{'filename':[star_file_ii]},
                                                },
                                            'normalization':{
                                                'type':'LuminosityNormalization',
                                                'SpecificLuminosityNormalization':{
                                                        'wavelength':norm_wl_ii,
                                                        'unitStyle':'neutralmonluminosity',
                                                        'specificLuminosity':norm_value_ii}
                                                    },
                                            'wavelengthBiasDistribution':self.__wavelengthBiasDistribution_options(ii,iteration0=iteration0)
                                            }
                                        }
            else:
                raise RuntimeError("This message should have not appeared!")
        
        
            translated_stellarSource = translate_dictionary(starSource_properties,0)
            source_ii  = translated_stellarSource[0]
            source_key = translated_stellarSource[1]
            number_ii = str(n_skirt_sources_count).zfill(3)
            #print("Stars: ", (source_key+" {}").format(number_ii))
            Sources['SourceSystem']['sources'][(source_key+" {}").format(number_ii)] = source_ii
    
        
        # GAS SOURCES
        #(skipped if iteration0 is True)
        
        if iteration0 == False:
            gasSource_properties = None
            
            all_gas_files       = get_from_folder(self._gas_sources_folder)
            gas_norm_wavelength = get_norm_wavelength(self._gas_sources_folder, 'nm')
            gas_norm_value      = get_norm(self._gas_sources_folder, 'erg/s')
            for ii in range(0,self._gas_zones):
                n_skirt_sources_count += 1
                
                gas_file_ii   = all_gas_files[ii]
                norm_wl_ii    = [gas_norm_wavelength[ii]]
                norm_value_ii = [gas_norm_value[ii]]
                
                if self._gas_geometry[ii][0] == 'shell':
                    Rin  = str(self._gas_geometry[ii][1])+self._geometryUnits #e.g.: '3.0 pc'
                    Rout = str(self._gas_geometry[ii][2])+self._geometryUnits
                    
                    gasSource_properties = {'GeometricSource':{
                                                'velocityMagnitude':'0 km/s', 
                                                'sourceWeight':"1", 
                                                'wavelengthBias': self._distributionBias,
                                                
                                                'geometry':{
                                                    'type':'Geometry',
                                                    'ShellGeometry':{
                                                        'minRadius': [Rin],
                                                        'maxRadius': [Rout], 
                                                        'exponent': ['0']
                                                        }
                                                    },
                                                'sed':{
                                                    'type':'SED',
                                                    'FileSED':{'filename':[gas_file_ii]}
                                                    },
                                                'normalization':{
                                                    'type':'LuminosityNormalization',
                                                    'SpecificLuminosityNormalization':{
                                                            'wavelength':norm_wl_ii,
                                                            'unitStyle':'neutralmonluminosity',
                                                            'specificLuminosity':norm_value_ii}
                                                    },
                                                'wavelengthBiasDistribution':self.__wavelengthBiasDistribution_options(ii,are_stars=False,iteration0=iteration0)
                                                }
                                            }
                    
                elif self._gas_geometry[ii][0] == 'ring':
                    R = str(self._gas_geometry[ii][1])+self._geometryUnits
                    w = str(self._gas_geometry[ii][2])+self._geometryUnits
                    h = str(self._gas_geometry[ii][3])+self._geometryUnits
                    
                    gasSource_properties = {'GeometricSource':{
                                                'velocityMagnitude':'0 km/s', 
                                                'sourceWeight':"1", 
                                                'wavelengthBias': self._distributionBias,
                                                
                                                'geometry':{
                                                    'type':'Geometry',
                                                    'RingGeometry':{
                                                        'ringRadius': [R], 
                            				            'width': [w], 
                            				            'height': [h]
                                                        }
                                                    },
                                                'sed':{
                                                    'type':'SED',
                                                    'FileSED':{'filename':[gas_file_ii]}
                                                    },
                                                'normalization':{
                                                    'type':'LuminosityNormalization',
                                                    'SpecificLuminosityNormalization':{
                                                            'wavelength':norm_wl_ii,
                                                            'unitStyle':'neutralmonluminosity',
                                                            'specificLuminosity':norm_value_ii}
                                                    },
                                                'wavelengthBiasDistribution':self.__wavelengthBiasDistribution_options(ii,are_stars=False,iteration0=iteration0)
                                                }
                                            }
                else:
                    raise RuntimeError("This message should not have appeared!")
            
            
                translated_gasSource = translate_dictionary(gasSource_properties,0)
                source_ii  = translated_gasSource[0]
                source_key = translated_gasSource[1]
                #number_ii = str(ii+1).zfill(3)
                number_ii = str(n_skirt_sources_count).zfill(3)
                #print("Gas: ", (source_key+" {}").format(number_ii))
                Sources['SourceSystem']['sources'][(source_key+" {}").format(number_ii)] = source_ii
        
        #returns
        self.Sources = Sources
        
    def __createMedia(self,iteration0):
        # =============================================================================
        # MEDIUM SYSTEM
        # =============================================================================
                
        Media = dict()
        
        # Basic properties 
        Media['MediumSystem'] = { 'photonPacketOptions':{'type':'PhotonPacketOptions',
                                                          'PhotonPacketOptions':{
                                                              'explicitAbsorption':'false',
                                                              'forceScattering':'true',
                                                              'minWeightReduction':1e4,
                                                              'minScattEvents':'0',
                                                              'pathLengthBias':'0.5'}
                                                         },
                                  'radiationFieldOptions':{'type':'RadiationFieldOptions',
                                                           'RadiationFieldOptions':{
                                                               'storeRadiationField': 'true',
                                                               'radiationFieldWLG':{'type':'DisjointWavelengthGrid',
                                                                    'LogWavelengthGrid':{
                                                                        'minWavelength': str(self._wavelength_min) + ' nm',
                                                                        'maxWavelength': str(self._wavelength_max) + ' nm',
                                                                        'numWavelengths': str(self._wavelength_res)
                                                                        }
                                                                    }
                                                                }
                                                          },
                                   #'dynamicStateOptions':{'type':'DynamicStateOptions',
                                   #                       'DynamicStateOptions':{'hasDynamicState':False, #If true, uncomment lines after 'recipes'
                                   #                                              'minIterations':'1',
                                   #                                              'maxIterations':'10',
                                   #                                              'iterationPacketsMultiplier':'1'
                                   #                                              }
                                   #                       }
                                }
        
        Media['MediumSystem']['media'] = {'type':'Medium'}
        
        # GAS CONTINUUM MEDIA (several mediums sharing the same geometry)
        
        gasMedia_properties = None
        #Define files needed for media if iteration0 is true or not
        Media_files = None
        Media_norm  = None
        if iteration0:
            Media_files = [self._nullMaterialFile for i in range(0,self._gas_zones)]
            Media_norm  = self._nullMass * np.ones(self._gas_zones)
        else:
            Media_files = get_from_folder(self._gas_opacity_folder)
            Media_norm  = get_mass_norm(self._gas_opacity_folder,'Msun') #It will check the third line of 'MeanFileDustMix', check if units are correct in your file
        
        for ii in range(0,self._gas_zones): 
        
            if self._gas_geometry[ii][0] == 'shell':
                Rin  = str(self._gas_geometry[ii][1])+self._geometryUnits #e.g.: '3.0 pc'
                Rout = str(self._gas_geometry[ii][2])+self._geometryUnits
                
                gasMedia_properties = {'GeometricMedium':{ #This line indicates the type of medium it is. Inside this dictionary you have ALL parameters needed for that medium
                                            'velocityMagnitude':'0 km/s', 
                                            'magneticFieldStrength':'0 uG',
                                            
                                            'geometry':{
                                                'type':'Geometry',
                                                'ShellGeometry':{
                                                    'minRadius': [Rin],  
                                                    'maxRadius': [Rout], 
                                                    'exponent': ['0']
                                                    }
                                                },
                                            'materialMix':{
                                                'type':'MaterialMix',
                                                'MeanFileDustMix':{'filename':[Media_files[ii]]}
                                                },
                                            'normalization':{
                                                'type':'MaterialNormalization',
                                                'MassMaterialNormalization':{
                                                        'mass':[Media_norm[ii]]
                                                }
                                            }
                                        }
                                    }
            elif self._gas_geometry[ii][0] == 'ring':
                R = str(self._gas_geometry[ii][1])+self._geometryUnits
                w = str(self._gas_geometry[ii][2])+self._geometryUnits
                h = str(self._gas_geometry[ii][3])+self._geometryUnits
                
                gasMedia_properties = {'GeometricMedium':{
                                            'velocityMagnitude':'0 km/s', 
                                            'magneticFieldStrength':'0 uG',
                                            
                                            'geometry':{
                                                'type':'Geometry',
                                                'RingGeometry':{
                                                    'ringRadius': [R],
                                                    'width': [w],
                                                    'height': [h]
                                                        }
                                                },
                                            'materialMix':{
                                                'type':'MaterialMix',
                                                'MeanFileDustMix':{'filename':[Media_files[ii]]}
                                                },
                                            'normalization':{
                                                'type':'MaterialNormalization',
                                                'MassMaterialNormalization':{
                                                        'mass':[Media_norm[ii]]
                                                }
                                            }
                                        }
                                    }
            else:
                raise RuntimeError("This message should really have never appeared")
        
        
            translated_gasMedia = translate_dictionary(gasMedia_properties,0)
            medium_ii  = translated_gasMedia[0]
            medium_key = translated_gasMedia[1]
            number_ii = str(ii+1).zfill(3)
            Media['MediumSystem']['media'][(medium_key+" {}").format(number_ii)] = medium_ii 
           

        Media['MediumSystem']['samplingOptions'] = {'type':'SamplingOptions',
                                                    'SamplingOptions':{
                                                        'numDensitySamples':100,
                                                        'numPropertySamples':1,
                                                        'aggregateVelocity':'Average'
                                                        }
                                                    }

        Media['MediumSystem']['grid'] = {'type':'SpatialGrid',
                           'Cylinder2DSpatialGrid': {'maxRadius':self._border_r, 
                                                     'minZ':'-'+self._border_z,
                                                     'maxZ':self._border_z,
                                                     'meshRadial':{'type':'Mesh',
                                                                   'LinMesh':{'numBins':self._resolution_r}},
                                                     'meshZ':{'type':'MoveableMesh',
                                                              'LinMesh':{'numBins':self._resolution_z}}
                                                     }
                           }       
        
        #returns
        self.Media = Media
            
    def __createInstruments(self):
        # =============================================================================
        # INSTRUMENTS
        # =============================================================================
        Instruments = dict()
         
        Instruments['InstrumentSystem'] = {'defaultWavelengthGrid':{
                'type':'WavelengthGrid',
                'LogWavelengthGrid':{
                    'minWavelength':str(self._wavelength_min)+' nm',
                    'maxWavelength':str(self._wavelength_max)+' nm',
                    'numWavelengths':str(self._wavelength_res)} 
                }
            }
        
        Instruments['InstrumentSystem']['instruments']={
                                'type':'Instrument',
                                 'SEDInstrument 1':{
                                         'instrumentName':'ised', 
                                         'distance':self._fluxProbeDistance,
                                         'inclination':'0 deg',
                                         'azimuth':'0',
                                         'roll':'0 deg',
                                         'radius':'0 pc',
                                         'recordComponents':'true', 
                                         'numScatteringLevels': '0', 
                                         'recordPolarization':'false', 
                                         'recordStatistics':'false',
                                         'wavelengthGrid':{'type':'WavelengthGrid',
                                                           'LogWavelengthGrid':{
                                                                   'minWavelength':str(self._wavelength_min)+' nm',
                                                                   'maxWavelength':str(self._wavelength_max)+' nm',
                                                                   'numWavelengths':str(self._wavelength_res)
                                                                               }
                                                   }
                                     },
                                 }
        #returns
        self.Instruments = Instruments
    
    def __createProbes(self):
        # =============================================================================
        # PROBES
        # =============================================================================
        
        Probes = dict()
        
        Probes['ProbeSystem'] = {
            'probes':{'type':'Probe',    
                'LuminosityProbe 1':{
                'probeName':'source_lum',
                    'wavelengthGrid':{'type':'WavelengthGrid',
                     'LogWavelengthGrid':{
                         'minWavelength':str(self._wavelength_min)+' nm',
                         'maxWavelength':str(self._wavelength_max)+' nm',
                         'numWavelengths':str(self._wavelength_res)
                                          }
                     }
                },
                'ConvergenceInfoProbe 1':{'probeName':"conv",
                                                 'wavelength':str(self._wavelength_norm)+' nm',
                                                 'probeAfter':'Setup'
                                                 },
                'RadiationFieldProbe 1':{
                    'probeName':"nuJnu",
                    'writeWavelengthGrid': "false",
                    'form':{
                        'type':'Form',
                        'AtPositionsForm':{
                            'filename': self.__Positions,
                            'useColumns':"",
                            }
                        },
                },
                'InstrumentWavelengthGridProbe 1':{'probeName':"grid"},
                'RadiationFieldWavelengthGridProbe 1':{'probeName':"radGrid"}
            }
            }
        #returns
        self.Probes = Probes
