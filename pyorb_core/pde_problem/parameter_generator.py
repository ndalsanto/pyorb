#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 14:21:47 2019

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np
import random

class Parameter_generator:
    """a class for handling the generation of parameters."""

    def __init__( self, _num_parameters):
        self.M_num_parameters = _num_parameters
        
    M_num_parameters = 0


class Random_parameter_generator( Parameter_generator ):
    """a class for handling the random generation of parameters."""

    def __init__( self, _num_parameters ):
        
        Parameter_generator.__init__( self, _num_parameters )

        return

    def generate_parameter( self, _mu_min, _mu_max ):
        
        new_parameter = np.zeros( self.M_num_parameters )
        
        for iP in range( new_parameter.shape[0] ):
            pRandom = float( random.randint(0,10000) ) / 10000.0
            new_parameter[iP] = _mu_min[iP] + pRandom * ( _mu_max[iP] - _mu_min[iP] )

        return new_parameter

class Tensor_parameter_generator( Parameter_generator ):
    """a class for handling the generation of tensor of parameters."""

    def __init__( self, _num_parameters, _parameter_steps ):
        
        Parameter_generator.__init__( self, _num_parameters )
        self.M_counter = np.zeros( self.M_num_parameters )
        self.M_parameter_steps = _parameter_steps
        
        return

    def generate_parameter( self, _mu_min, _mu_max ):
        
        new_parameter = np.zeros( self.M_num_parameters )
        
        for iP in range( self.M_num_parameters ):
            new_parameter[iP] = _mu_min[iP] + ( _mu_max[iP] - _mu_min[iP] ) \
                                            * self.M_counter[iP].astype(float) / ( self.M_parameter_steps[iP].astype(float) - 1)

        self.update_next_parameter( )
        return new_parameter

    def update_next_parameter( self ):
        
        current_update = 0
        
        while current_update < self.M_num_parameters and self.M_counter[ current_update ] >= self.M_parameter_steps[ current_update ]-1:
            self.M_counter[ current_update ] = 0
            current_update += 1
            
        if current_update < self.M_num_parameters:
            self.M_counter[ current_update ] += 1

        return
        
    def reset_counter( self ):
        self.M_counter = np.zeros( self.M_num_parameters )
        
        
    
    M_current_update = 0
    M_counter = np.zeros( 0 )
    M_parameter_steps = np.zeros( 0 )

