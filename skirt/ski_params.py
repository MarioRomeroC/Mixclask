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
    def __init__(self,gas_path,star_path,output_positions,wavelength_dict):
        self.__initParams()
        self.__initWavelengths(wavelength_dict)
        self.__initFolders()
        self.__initDetails()
        
        self.__parseData(gas_path,'Gas')
        self.__parseData(star_path,'Star')
        
        self.__deduceLimits()
        self.__checkSkirtIssues()
        self.__MeanIntensity_Positions = output_positions
    
    ### INITIALIZATION
    
    def __initParams(self):
        #Geometry GAS params
        self._gas_geometry  = None
        #'shell' geometry
        self._gas_minRadius = []
        self._gas_maxRadius = []
        #'ring' geometry
        self._gas_ringRadius = []
        self._gas_ringWidth  = []
        self._gas_ringHeight = []
        #number of gas zones
        self._gas_zones = None
        
        #Geometry STAR params
        self._star_geometry = None
        self._star_files = []
        #'shell' geometry
        self._star_minRadius = []
        self._star_maxRadius = []
        #'ring' geometry
        self._star_ringRadius = []
        self._star_ringWidth  = []
        self._star_ringHeight = []
        #'point' geometry
        self._star_x = []
        self._star_y = []
        self._star_z = []
        #number of stellar zones
        self._star_zones = None
    
    def __initFolders(self):   
        #folder locations
        # DO NOT CHANGE THESE PARAMETERS, THEY ARE HARDCODED BECAUSE CLOUDY MOVES THE FILES THERE
        self._gas_opacity_folder  = "input_data/gas_props"
        self._gas_sources_folder  = "input_data/gas_sources"
        # Except next one, you can locate this wherever you want
        self._star_sources_folder = "input_data/star_sources"
    
    def __initWavelengths(self,wavelength_dict):
        self._wavelength_max  = wavelength_dict['maxWavelength']
        self._wavelength_min  = wavelength_dict['minWavelength']
        self._wavelength_norm = wavelength_dict['normalization']
        self._wavelength_res  = wavelength_dict['resolution']
    
    def __initDetails(self):
        #This function contains parameters used for debugging
        self._photonPackets  = '1e8'
        self._montecarloSeed = '0'
        self._fluxProbeDistance = '1 Mpc'
        self._geometryUnits = ' pc' #The space IS intended
        self._resolution_r = 200
        self._resolution_z = 200
        #Skirt cries if you use 'ring' geometry and place one at R = 0pc
        self._skirt_ringRadius_correction = 0.01 #pc
        #Iteration0 options
        self._nullMaterialFile = 'params/NullMaterial.stab' #Skirt will use this file as 'material' in the iteration0 (should be a material without absorption and scattering)
        self._nullMass = 1.0e-10 #Msun. Normalization of 'NullMaterial'. The ideal value is 0.0, but I selected a small number because skirt will not run instead

    
    def __deduceLimits(self):
        #limits are decided with gas zones.
        if self._gas_geometry == 'shell':
            self._border_r = self._gas_maxRadius[-1] #This should give 'X pc' or similar, depending of 'self._geometry_units'
            self._border_z = self._border_r
        elif self._gas_geometry == 'ring':
            #Not that straightforward
            self._border_r  = float(self._gas_ringRadius[-1].replace(self._geometryUnits,''))
            self._border_r += float(self._gas_ringWidth[-1].replace(self._geometryUnits,''))
            self._border_r  = str(self._border_r)+self._geometryUnits
            #Especially z-axis limits.
            factor = 5.0 #Because potato
            self._border_z  = float(self._gas_ringHeight[-1].replace(self._geometryUnits,''))
            self._border_z *= factor
            self._border_z  = str(self._border_z)+self._geometryUnits
    
    def __checkSkirtIssues(self):
        #(1) If ring geometry is set and some of 'ringRadius' is 0.0
        #   Skirt crashes due they do not accept rings with R=0
        radius_correction = '0.01 pc'
        if self._gas_geometry == 'ring':# and any(R == '0.0 pc' for R in self._gas_ringRadius):
            for r in range(0,len(self._gas_ringRadius)):
                if self._gas_ringRadius[r] == '0.0 pc':
                    self._gas_ringRadius[r] = radius_correction
        #Same issue with stars
        if self._star_geometry == 'ring':
            for r in range(0,len(self._star_ringRadius)):
                if self._star_ringRadius[r] == '0.0 pc':
                    self._star_ringRadius[r] = radius_correction
    
    ### PARSING DATA
    
    def __parseData(self,file_path,Type):
        #This goes a little different than 'parseData' in CloudyClass.py
        #Type has two options: 'Gas' or 'Star'
        file = open(file_path,'r')
        
        #Prepare some arrays first to avoid too many if evaluations
        file_geometry = None
        
        file_minRadius = []
        file_maxRadius = []
        
        file_ringRadius = []
        file_ringWidth  = []
        file_ringHeight = []
        
        file_x = []
        file_y = []
        file_z = []
        
        #Read the file
        geometry_units = self._geometryUnits
        n_rows = 0
        starting_col = 1 if Type == 'Star' else 3 #For gas files, first columns are the sed file, the mass and nH
        while True:
            line = file.readline()
            if not line: break #EoF
            line_data = line.split()
            
            if line_data[0] == '#' and line_data[1].lower() == 'geometry':
                #parse options
                #line for parsing is '# geometry : shell', therefore
                file_geometry = line_data[3].lower()
            elif line_data[0] == '#': 
                continue
            else:
                n_rows += 1
                #Get relevant data
                #Column 1 is always the sed file for stars:
                if Type == 'Star':
                    self._star_files.append(line_data[0])
                #Now the variable columns
                col = starting_col
                ## Geometry related columns
                if file_geometry == 'shell':
                    file_minRadius.append(line_data[col]+geometry_units)
                    file_maxRadius.append(line_data[col+1]+geometry_units)
                    col += 2
                elif file_geometry == 'ring':
                    file_ringRadius.append(line_data[col]+geometry_units)
                    file_ringWidth.append(line_data[col+1]+geometry_units)
                    file_ringHeight.append(line_data[col+2]+geometry_units)
                    col += 3
                elif file_geometry == 'point':
                    file_x.append(line_data[col]+geometry_units)
                    file_y.append(line_data[col+1]+geometry_units)
                    file_z.append(line_data[col+2]+geometry_units)
                else:
                    raise RuntimeError("Geometry not implemented!")
                    
        #And the final param
        file_zones = n_rows
        
        file.close()
        
        #Now tell me if data is 'gas' or 'star'
        if Type == 'Gas':
            self._gas_geometry = file_geometry
            self._gas_zones = file_zones
            if self._gas_geometry == 'shell':
                self._gas_minRadius = file_minRadius
                self._gas_maxRadius = file_maxRadius
            elif self._gas_geometry == 'ring':
                self._gas_ringRadius = file_ringRadius
                self._gas_ringWidth  = file_ringWidth
                self._gas_ringHeight = file_ringHeight
            else:
                raise RuntimeError("Gas can never be a point source!")
        elif Type == 'Star':
            self._star_geometry = file_geometry
            self._star_zones = file_zones
            if self._star_geometry == 'shell':
                self._star_minRadius = file_minRadius
                self._star_maxRadius = file_maxRadius
            elif self._star_geometry == 'ring':
                self._star_ringRadius = file_ringRadius
                self._star_ringWidth  = file_ringWidth
                self._star_ringHeight = file_ringHeight
            elif self._star_geometry == 'point':
                self._star_x = file_x
                self._star_y = file_y
                self._star_z = file_z
            else:
                raise RuntimeError("This error message should have not appeared!")
        
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
                        'numPackets':self._photonPackets,
                    'random':{
                        'type':'Random',
                        'Random':{'seed':self._montecarloSeed}
                                },        
                    'units':{
                        'type':'Units',
                        'ExtragalacticUnits':{'fluxOutputStyle': 'Neutral'}
                            },
                    'cosmology':{
                        'type':'Cosmology',
                        'LocalUniverseCosmology':{}
                                }
                    }
                }
        #Return
        self.Basics = Basics
    
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
        if self._star_geometry == 'shell':
            starSource_properties = {'GeometricSource':{ #This line indicates the type of source it is. Inside this dictionary you have ALL parameters needed for that source
                		                'velocityMagnitude':'0 km/s', 
                		                'sourceWeight':"1", 
                		                'wavelengthBias':"0.5",
                		                
                		                'geometry':{
                		                    'type':'Geometry',
                		                    'ShellGeometry':{
                    				            'minRadius': self._star_minRadius, #All lists (here labelled as a variable) here must have length equal to 'n_gasSources'. 
                    				            'maxRadius': self._star_maxRadius, 
                    				            'exponent': ['0']*self._star_zones
                    		                        }
                		                    },
                		                'sed':{
                		                    'type':'SED',
                		                    'FileSED':{'filename':get_from_folder(self._star_sources_folder, self._star_files)}
                		                    },
                		                'normalization':{
                		                    'type':'LuminosityNormalization',
                		                    'SpecificLuminosityNormalization':{
                		                            'wavelength':get_norm_wavelength(self._star_sources_folder, 'nm', self._star_files), #It will check the third line of 'fileSED', check if units are correct in your file
                		                            'unitStyle':'neutralmonluminosity',
                		                            'specificLuminosity':get_norm(self._star_sources_folder, 'erg/s', self._star_files)} #It will check the fourth line of 'fileSED',check if units are correct in your file
                		                    }
                		                }
                		            }
        elif self._star_geometry == 'ring':
            starSource_properties = {'GeometricSource':{ #This line indicates the type of source it is. Inside this dictionary you have ALL parameters needed for that source
                		                'velocityMagnitude':'0 km/s', 
                		                'sourceWeight':"1", 
                		                'wavelengthBias':"0.5",
                		                
                		                'geometry':{
                		                    'type':'Geometry',
                		                    'RingGeometry':{
                    				            'ringRadius': self._star_ringRadius, #All lists (here labelled as a variable) here must have length equal to 'n_gasSources'. 
                    				            'width': self._star_ringWidth, 
                    				            'height': self._star_ringHeight
                		                        }
                		                    },
                		                'sed':{
                		                    'type':'SED',
                		                    'FileSED':{'filename':get_from_folder(self._star_sources_folder, self._star_files)}
                		                    },
                		                'normalization':{
                		                    'type':'LuminosityNormalization',
                		                    'SpecificLuminosityNormalization':{
                		                            'wavelength':get_norm_wavelength(self._star_sources_folder, 'nm', self._star_files), #It will check the third line of 'fileSED', check if units are correct in your file
                		                            'unitStyle':'neutralmonluminosity',
                		                            'specificLuminosity':get_norm(self._star_sources_folder, 'erg/s', self._star_files)} #It will check the fourth line of 'fileSED',check if units are correct in your file
                		                    }
                		                }
                		            }
        elif self._star_geometry == 'point':
            starSource_properties = {'PointSource':{
                                        'positionX': self._star_x,
                                        'positionY': self._star_y,
                                        'positionZ': self._star_x,
                                        'velocityX': ['0 km/s']*self._star_zones,
                                        'velocityY': ['0 km/s']*self._star_zones,
                                        'velocityZ': ['0 km/s']*self._star_zones,
                                        'sourceWeight': '1',
                                        'wavelengthBias':'0.5',
                                        'angularDistribution': {
                                            'type':'AngularDistribution',
                                            'IsotropicAngularDistribution':None #This is to copy '{}' in the dictionary
                                            },
                                        'polarizationProfile': {
                                            'type': 'PolarizationProfile',
                                            'NoPolarizationProfile':None
                                            },
                                            
                                        'sed':{ 'type':'SED',
                                            'FileSED':{'filename':get_from_folder(self._star_sources_folder, self._star_files)}},
                                        'normalization':{
                                            'type':'LuminosityNormalization',
                                            'SpecificLuminosityNormalization':{
                                                    'wavelength':get_norm_wavelength(self._star_sources_folder, 'nm', self._star_files),
                                                    'unitStyle':'neutralmonluminosity',
                                                    'specificLuminosity':get_norm(self._star_sources_folder, 'erg/s', self._star_files)}
                                            }
                                        }
                                    }
        else:
            raise RuntimeError("This message should have not appeared!")
        
        for ii in range(self._star_zones):
            n_skirt_sources_count += 1
            translated_stellarSource = translate_dictionary(starSource_properties,ii)
            source_ii  = translated_stellarSource[0]
            source_key = translated_stellarSource[1]
            number_ii = str(n_skirt_sources_count).zfill(3)
            #print("Stars: ", (source_key+" {}").format(number_ii))
            Sources['SourceSystem']['sources'][(source_key+" {}").format(number_ii)] = source_ii
    
        
        # GAS SOURCES
        #(skipped if iteration0 is True)
        
        if iteration0 == False:
            gasSource_properties = None
            #Mario: 0.0 pc is not allowed by skirt
            if self._gas_geometry == 'shell':
                gasSource_properties = {'GeometricSource':{
                                            'velocityMagnitude':'0 km/s', 
                                            'sourceWeight':"1", 
                                            'wavelengthBias':"0.5",
                                            
                                            'geometry':{
                                                'type':'Geometry',
                                                'ShellGeometry':{
                                                    'minRadius': self._gas_minRadius,
                                                    'maxRadius': self._gas_maxRadius, 
                                                    'exponent': ['0']*self._gas_zones
                                                    }
                                                },
                                            'sed':{
                                                'type':'SED',
                                                'FileSED':{'filename':get_from_folder(self._gas_sources_folder)}
                                                },
                                            'normalization':{
                                                'type':'LuminosityNormalization',
                                                'SpecificLuminosityNormalization':{
                                                        'wavelength':get_norm_wavelength(self._gas_sources_folder, 'nm'),
                                                        'unitStyle':'neutralmonluminosity',
                                                        'specificLuminosity':get_norm(self._gas_sources_folder, 'erg/s')}
                                                }
                                            }
                                        }
            elif self._gas_geometry == 'ring':
                gasSource_properties = {'GeometricSource':{
                                            'velocityMagnitude':'0 km/s', 
                                            'sourceWeight':"1", 
                                            'wavelengthBias':"0.5",
                                            
                                            'geometry':{
                                                'type':'Geometry',
                                                'RingGeometry':{
                                                    'ringRadius': self._gas_ringRadius, #All lists (here labelled as a variable) here must have length equal to 'n_gasSources'. 
                        				            'width': self._gas_ringWidth, 
                        				            'height': self._gas_ringHeight
                                                    }
                                                },
                                            'sed':{
                                                'type':'SED',
                                                'FileSED':{'filename':get_from_folder(self._gas_sources_folder)}
                                                },
                                            'normalization':{
                                                'type':'LuminosityNormalization',
                                                'SpecificLuminosityNormalization':{
                                                        'wavelength':get_norm_wavelength(self._gas_sources_folder, 'nm'),
                                                        'unitStyle':'neutralmonluminosity',
                                                        'specificLuminosity':get_norm(self._gas_sources_folder, 'erg/s')}
                                                }
                                            }
                                        }
            else:
                raise RuntimeError("This message should not have appeared!")
            
            for ii in range(self._gas_zones):
                n_skirt_sources_count += 1
                translated_gasSource = translate_dictionary(gasSource_properties,ii)
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
            
            #Define files needed for media if iteration0 is true or not
            Media_files = None
            Media_norm  = None
            if iteration0:
                Media_files = [self._nullMaterialFile for i in range(0,len(self._gas_zones))]
                Media_norm  = self._nullMass * np.ones(self._gas_zones)
            else:
                Media_files = get_from_folder(self._gas_opacity_folder)
                Media_norm  = get_mass_norm(self._gas_opacity_folder,'Msun') #It will check the third line of 'MeanFileDustMix', check if units are correct in your file
            
            # Basic properties 
            Media['MediumSystem'] = { 'numDensitySamples':'100',
                                       'photonPacketOptions':{'type':'PhotonPacketOptions',
                                                              'PhotonPacketOptions':{
                                                                  'forceScattering':True,
                                                                  'minWeightReduction':1e4,
                                                                  'minScattEvents':'0',
                                                                  'pathLengthBias':'0.5'}},
                   'extinctionOnlyOptions':{'type':'ExtinctionOnlyOptions',
                                            'ExtinctionOnlyOptions':{
                                                'storeRadiationField':'true',
                                                'radiationFieldWLG':{'type':'DisjointWavelengthGrid',
                                                                 'LogWavelengthGrid':{
                                                                     'minWavelength':str(self._wavelength_min)+' nm',
                                                                     'maxWavelength':str(self._wavelength_max)+' nm',
                                                                     'numWavelengths': str(self._wavelength_res)
                                                                                     }
                                                                 }}
                                            },
                   'dynamicStateOptions':{'type':'DynamicStateOptions', 
                                          'DynamicStateOptions':{'hasDynamicState':False, #If true, uncomment lines after 'recipes'
                                                                 'minIterations':'1',
                                                                 'maxIterations':'10',
                                                                 'iterationPacketsMultiplier':'1'
                                                                 }
                                          }
                                        }
            
            Media['MediumSystem']['media'] = {'type':'Medium'}
            
            # GAS CONTINUUM MEDIA (several mediums sharing the same geometry)
            
            gasMedia_properties = None
            
            if self._gas_geometry == 'shell':
                gasMedia_properties = {'GeometricMedium':{ #This line indicates the type of medium it is. Inside this dictionary you have ALL parameters needed for that medium
                                            'velocityMagnitude':'0 km/s', 
                                            'magneticFieldStrength':'0 uG',
                                            
                                            'geometry':{
                                                'type':'Geometry',
                                                'ShellGeometry':{
                                                    'minRadius': self._gas_minRadius, #All lists (here labelled as a variable) here must have length equal to 'n_gasSources'. 
                                                    'maxRadius': self._gas_maxRadius, 
                                                    'exponent': ['0']*self._gas_zones
                                                    }
                                                },
                                            'materialMix':{
                                                'type':'MaterialMix',
                                                'MeanFileDustMix':{'filename':Media_files}
                                                },
                                            'normalization':{
                                                'type':'MaterialNormalization',
                                                'MassMaterialNormalization':{
                                                        'mass':Media_norm
                                                }
                                            }
                                        }
                                    }
            elif self._gas_geometry == 'ring':
                gasMedia_properties = {'GeometricMedium':{
                                            'velocityMagnitude':'0 km/s', 
                                            'magneticFieldStrength':'0 uG',
                                            
                                            'geometry':{
                                                'type':'Geometry',
                                                'RingGeometry':{
                                                    'ringRadius': self._gas_ringRadius,
                                                    'width': self._gas_ringWidth,
                                                    'height': self._gas_ringHeight
                                                        }
                                                },
                                            'materialMix':{
                                                'type':'MaterialMix',
                                                'MeanFileDustMix':{'filename':Media_files}
                                                },
                                            'normalization':{
                                                'type':'MaterialNormalization',
                                                'MassMaterialNormalization':{
                                                        'mass':Media_norm
                                                }
                                            }
                                        }
                                    }
            else:
                raise RuntimeError("This message should really have never appeared")
            
            for ii in range(self._gas_zones): 
                translated_gasMedia = translate_dictionary(gasMedia_properties,ii)
                medium_ii  = translated_gasMedia[0]
                medium_key = translated_gasMedia[1]
                number_ii = str(ii+1).zfill(3)
                Media['MediumSystem']['media'][(medium_key+" {}").format(number_ii)] = medium_ii 
               
                   
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
            'SpatialGridConvergenceProbe 1':{'probeName':"conv",
                                             'wavelength':str(self._wavelength_norm)+' nm'
                                             },
            'RadiationFieldAtPositionsProbe 1':{
                'probeName':"nuJnu",
                'filename':self.__MeanIntensity_Positions,
                'useColumns':"",
                'writeWavelengthGrid':"false"
                                            },
            'InstrumentWavelengthGridProbe 1':{'probeName':"grid"},    
            'RadiationFieldWavelengthGridProbe 1':{'probeName':"radGrid"}                        
                    }
            }
        #returns
        self.Probes = Probes