#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:58:02 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np

import fom_problem as fp

def tb_theta_a( _param, _q ):
    
    assert _q <= 3
    
    if( _q < 3 ):
        return _param[_q]
    else:
        return 1.0

def tb_theta_f( _param, _q ):
    return 1.0


class thermal_block_problem( fp.fom_problem ):

    def __init__( self ):
        fp.fom_problem.__init__( self )
        
        return
   
    def define_theta_functions( self ):
        
        self.M_theta_a = tb_theta_a
        self.M_theta_f = tb_theta_f
        
        return
     
    def solve_external_fom_problem( self, _param ):
        
        print( "Solving the FOM problem with parameter \n" )
        print( _param ) 
        
        converted_param = self.convert_parameter( _param )
             
        print( "Converted parameter is \n" )
        print( converted_param ) 
           
        sol = self.M_external_engine.solve_parameter( converted_param, self.M_fom_specifics )

        return sol
        
  
    
