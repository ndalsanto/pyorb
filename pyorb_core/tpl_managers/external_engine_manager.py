#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 13:59:45 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

from pyorb_core.tpl_managers import matlab_external_engine as matlab_ext
from pyorb_core.tpl_managers import cpp_external_engine as cpp_ext

class external_engine_manager( ):
    
    def __init__( self, _engine_type, _library_path ):
        
        self.M_engine_type  = _engine_type
        self.M_library_path = _library_path
        
        if _engine_type == 'matlab':
            self.M_external_engine = matlab_ext.matlab_external_engine( _engine_type, _library_path )
        elif _engine_type == 'cpp':
            self.M_external_engine = cpp_ext.cpp_external_engine( _engine_type, _library_path )
        
        return

    M_engine_type = ""
    M_library_path = ""
    M_external_engine = None
    
    def get_external_engine( self ):
        return self.M_external_engine

    def start_engine( self ):
        self.M_external_engine.start_engine( )
        return
    
    def quit_engine( self ):
        self.M_external_engine.quit_engine( )
        return




