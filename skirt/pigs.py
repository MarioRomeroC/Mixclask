#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  7 15:20:20 2021

@author: pablo
"""

import importlib
from glob import glob
from numpy import sort
from os import path, getcwd

class SkirtFile(object):
    """
    This class generates a .ski file for running SKIRT models with
    multiple built-in geometry configurations. 
    
    input: 
        - config file (name of the .py configuration file specifiying all the
                        parameters under request)
        - output name (Name or full path of the resultant .ski file)    
    """
    
    
    def __init__(self, params, output_path=None):
        # self.config_file = config_file_path                
        self.params = params        
        if output_path:        
            output = output_path+'.ski'
        else:
            output = path.join(getcwd(), 'skirt_config_file.ski')
            
        self.skifile = open(output, 'w+')
            
        self.skifile.write(
            '<?xml version="1.0" encoding="UTF-8"?>\n<!-- A SKIRT parameter file created by PIGS! -->\n'
            )                

        self.make_basics()
        self.make_sources()
        self.make_media()
        self.make_instruments()
        self.make_probes()
        
        self.skifile.write(             
        '    </MonteCarloSimulation>\n</skirt-simulation-hierarchy>'
                            )

        print('File saved at ', output)
    # def read_config(self):
    #     cfg = importlib.import_module(self.config_file)        
    #     return cfg
    
    def make_basics(self):        
        first_line = '<skirt-simulation-hierarchy type="MonteCarloSimulation" format="9" producer="SKIRT v9.0" time="---">\n'
        end_line = ''
        lines = self.write_dict(dictionary=self.params.Basics)        
        basics_code = first_line+lines[3:]+end_line
             
        basics_code = basics_code.replace(' />', '/>')
        basics_code = basics_code.replace(' >', '>')
        basics_code = basics_code.replace('</MonteCarloSimulation>', '')
        self.skifile.write(basics_code)
    
    def make_sources(self):
        first_line = '<sourceSystem type="SourceSystem">\n'
        end_line = '\n        </sourceSystem>'
        lines = self.write_dict(dictionary=self.params.Sources)        
        sources_code = first_line+lines[3:]+end_line
        
        all_keys = list(self.params.Sources['SourceSystem']['sources'].keys())
        for key_i in all_keys:
            if key_i != 'type':
                key_split = key_i.split()            
                sources_code = sources_code.replace(key_i, key_split[0])                   
                
        sources_code = sources_code.replace(' />', '/>')        
        sources_code = sources_code.replace(' >', '>')        
        self.skifile.write(sources_code)
        
                                    
    def make_media(self):
        first_line = '<mediumSystem type="MediumSystem">\n'
        end_line = '\n        </mediumSystem>'
        lines = self.write_dict(dictionary=self.params.Media)        
        media_code = first_line+lines[3:]+end_line
        
        all_keys = list(self.params.Media['MediumSystem']['media'].keys())
        for key_i in all_keys:
            if key_i != 'type':
                key_split = key_i.split()            
                media_code = media_code.replace(key_i, key_split[0])                   
        media_code = media_code.replace(' />', '/>')        
        media_code = media_code.replace(' >', '>')        
        self.skifile.write(media_code)
                        
    def make_instruments(self):
        
        first_line = '<instrumentSystem type="InstrumentSystem">\n'
        end_line = '\n        </instrumentSystem>'
        lines = self.write_dict(dictionary=self.params.Instruments)        
        instruments_code = first_line+lines[3:]+end_line
        
        all_keys = list(
            self.params.Instruments['InstrumentSystem']['instruments'].keys())
        for key_i in all_keys:
            if key_i != 'type':
                key_split = key_i.split()            
                instruments_code = instruments_code.replace(key_i, key_split[0])                   
        instruments_code = instruments_code.replace(' />', '/>') 
        instruments_code = instruments_code.replace('  >', '>') 
        self.skifile.write(instruments_code)
                    
    def make_probes(self):
        first_line = '<probeSystem type="ProbeSystem">\n'
        end_line = '\n        </probeSystem>'
        lines = self.write_dict(dictionary=self.params.Probes)        
        probes_code = first_line+lines[3:]+end_line
        
        all_keys = list(
            self.params.Probes['ProbeSystem']['probes'].keys())
        for key_i in all_keys:
            if key_i != 'type':
                key_split = key_i.split()            
                probes_code = probes_code.replace(key_i, key_split[0])                   
        probes_code = probes_code.replace(' />', '/>')        
        probes_code = probes_code.replace(' >', '>')        
        self.skifile.write(probes_code)
           
     
    def write_dict(self, dictionary, name=None, identation=''):        
        lines = ''
        has_child = False
        
        for key, value in dictionary.items():
            if isinstance(value, dict):                
                if not lines[-1:] == '>':                    
                    lines +='>\n    '
                identation += '    '
                lines += self.write_dict(value, key, identation=identation)
                has_child = True
            else:
                lines += '{}="{}" '.format(key, value)        
        if name:
            header = identation+'<{} '.format(name)
            if has_child:
                ending = '\n</{}>'.format(name)            
            else:
                ending = '/>'
        else:
            header = ''
            ending = ''
            
        lines = header+lines+ending
        return lines
           
        
# ...

if __name__ == '__main__':
    skirtfile = SkirtFile('ski_params', 
                          output_path='my_skirt_config_file.ski')