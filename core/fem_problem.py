#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 12:02:21 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

def default_theta_function( _param, _q ):

    print( "You are using the default theta function, please provide specific ones for your problem " )
    
    pass

class fem_problem( ):
    
    def __init__( self ):
        return
    
    def get_theta_a( self, _param, _q ):
        return self.M_theta_a( _param, _q )
    
    def get_theta_f( self, _param, _q ):
        return self.M_theta_f( _param, _q )
    
    M_theta_a = default_theta_function
    M_theta_f = default_theta_function   
    
    