#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:49:46 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np

class Parameter_handler:
    """a class for handling the parameters."""

    def __init__( self, _parameter_generator ):

        self.M_parameter_generator = _parameter_generator
        return

    def assign_parameters_bounds( self, _param_min, _param_max ):
        self.M_param_min = _param_min
        self.M_param_max = _param_max
        self.M_param     = np.zeros( _param_min.shape )
        self.M_num_parameters = _param_min.shape[0]

    def assign_parameters( self, _param ):
        
        assert self.M_num_parameters == _param.shape[0]
        self.M_param = _param

    def print_parameters( self ):
        print( "Numberof parameters : %d " % self.M_num_parameters )
        print( "The current parameter is: " )
        print( self.M_param )

    def substitute_parameter_generator( self, _parameter_generator ):

        self.M_parameter_generator = _parameter_generator
        return

    def generate_parameter( self ):
        assert self.M_num_parameters > 0
        self.M_param = self.M_parameter_generator.generate_parameter( self.M_param_min, self.M_param_max )
        return
        
    def get_parameter( self ):
        return self.M_param

    def get_parameter_vector( self ):
        return self.M_param

    def get_num_parameters( self ):
        return self.M_num_parameters

    def get_min_parameters( self ):
        return self.M_param_min

    def get_max_parameters( self ):
        return self.M_param_max

    def get_range_parameters( self ):
        return self.M_param_max - self.M_param_min

    def set_parameter_generator( self, _parameter_generator ):
        self.M_parameter_generator = _parameter_generator
        
    M_param_min = np.zeros( 0 )
    M_param_max = np.zeros( 0 )
    M_param     = np.zeros( 0 )
    M_num_parameters = 0
    M_parameter_generator = None

