#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 11:58:02 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""






def theta_a( _param, _q ):
    
    assert _q <= 3
    
    if( _q < 3 ):
        return _param[_q]
    else:
        return 1.0

def theta_f( _param, _q ):
    return 1.0

