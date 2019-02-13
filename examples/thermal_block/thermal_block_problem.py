#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:58:02 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np

import pyorb_core.pde_problem.fom_problem as fp

def tb_theta_a( _param, _q ):
    
    assert _q <= 3
    
    if( _q < 3 ):
        return _param[_q]
    else:
        return 1.0

def tb_theta_f( _param, _q ):
    return 1.0

def tb_full_theta_f( _param ):
    return np.array([1.0])

def tb_full_theta_a( _param ):
    
    diffusions = np.zeros( len( _param) + 1 )
    
    diffusions[0:len( _param)] = _param
    
    diffusions[len( _param)] = 1.0
    
    return diffusions


class thermal_block_problem( fp.fom_problem ):

    def __init__( self, _parameter_handler ):
        fp.fom_problem.__init__( self, _parameter_handler )
        
        return
   
    def define_theta_functions( self ):
        
        self.M_theta_a = tb_theta_a
        self.M_theta_f = tb_theta_f
        
        self.M_full_theta_a = tb_full_theta_a
        self.M_full_theta_f = tb_full_theta_f

        return
     

  
    
