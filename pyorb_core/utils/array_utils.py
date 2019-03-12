#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 24 10:23:20 2019

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np

def save_matrix( _matrix, _file_name ):
    
    output_file = open( _file_name, 'w+' )
    
    for iNs in range( _matrix.shape[0] ):
        for iQ in range( _matrix.shape[1] ):
            output_file.write( "%.10g" % _matrix[iNs, iQ] )

            if iQ < _matrix.shape[1] - 1:
                output_file.write( " " )
            else:
                output_file.write( "\n" )
    
    output_file.close( )

def save_vector( _matrix, _file_name ):
    
    output_file = open( _file_name, 'w+' )
    
    for iNs in range( _matrix.shape[0] ):
        output_file.write( "%.10g " % _matrix[iNs] )
    
    output_file.close( )
