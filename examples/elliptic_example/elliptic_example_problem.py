#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:58:02 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np

import pyorb_core.pde_problem.fom_problem as fp

def ell_ex_theta_a( _param, _q ):
    
    assert _q <= 1 and _q >=0
    
    if( _q == 0 ):
        return 1.0
    elif _q == 1:
        return _param[0]

def ell_ex_theta_f( _param, _q ):

    assert _q <= 1 and _q >=0
    
    if( _q == 0 ):
        return _param[1]
    elif _q == 1:
        return _param[2]


class elliptic_example_problem( fp.fom_problem ):

    def __init__( self, _parameter_handler ):
        fp.fom_problem.__init__( self, _parameter_handler )
        
        return
   
    def define_theta_functions( self ):
        
        self.M_theta_a = ell_ex_theta_a
        self.M_theta_f = ell_ex_theta_f
        
        return
     

    
