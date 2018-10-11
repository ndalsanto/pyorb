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

    def __init__( self ):
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

    def generate_parameter( self ):
        # generate numbers between 0 and 1
        assert self.M_num_parameters > 0
        
        for iP in range( self.M_num_parameters ):
            pRandom = float( random.randint(0,10000) ) / 10000.0
            self.M_param[iP] = self.M_param_min[iP] + pRandom * ( self.M_param_max[iP] - self.M_param_min[iP] )
        
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


    M_param_min = np.zeros( 0 )
    M_param_max = np.zeros( 0 )
    M_param     = np.zeros( 0 )
    M_num_parameters = 0
