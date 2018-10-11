#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:58:02 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import fem_problem as fp

def tb_theta_a( _param, _q ):
    
    assert _q <= 3
    
    if( _q < 3 ):
        return _param[_q]
    else:
        return 1.0

def tb_theta_f( _param, _q ):
    return 1.0



class thermal_block_problem( fp.fem_problem ):

    def __init__( self ):
        fp.fem_problem.__init__( self )
        
        self.M_theta_a = tb_theta_a
        self.M_theta_f = tb_theta_f
        
        return
   
    
