#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:17:29 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""


import sys
import matlab.engine
#sol = eng.solve_parameter( matlab.double([param]), n_elements_x )


class external_engine_manager( ):
    
    def __init__( self, _engine_type, _library_path ):
        
        self.M_engine_type  = _engine_type
        self.M_library_path = _library_path
        
        return

    M_engine_type = ""
    M_library_path = ""
    M_engine = 0
    
    def get_engine( self ):
        
        return self.M_engine
    
    def start_engine( self ):
        
        engine_switcher = {
                    'matlab': self.start_matlab_engine # , \
                }
        
        # retrieving the function activating the correct engine and giving a lambda as default value
        engine_activation = engine_switcher.get( self.M_engine_type, lambda: "Invalid engine activation" )
        
        engine_activation( )
        
        return
    
    def quit_engine( self ):
        
        engine_switcher = {
                    'matlab': self.quit_matlab_engine # , \
                }
        
        # retrieving the function activating the correct engine and giving a lambda as default value
        engine_deactivation = engine_switcher.get( self.M_engine_type, lambda: "Invalid engine deactivation" )
        
        engine_deactivation( )
        
        return

    
    def start_matlab_engine( self ):
        
        self.M_engine = matlab.engine.start_matlab( )
        self.M_engine.addpath( self.M_engine.genpath( self.M_library_path ) )

        print( 'Successfully started matlab engine and corresponding FOM library' )

        return

    def quit_matlab_engine( self ):
        
        self.M_engine.quit( )
        
        print( 'Successfully quitted matlab engine' )

        return




