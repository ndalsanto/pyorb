#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 18:43:23 2018

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


class nonaffine_diffusion_problem( fp.fom_problem ):

    def __init__( self, _parameter_handler ):
        fp.fom_problem.__init__( self, _parameter_handler )
        
        return
   
    def define_theta_functions( self ):
        
        self.M_theta_a = tb_theta_a
        self.M_theta_f = tb_theta_f
        
        return
     

  
    
