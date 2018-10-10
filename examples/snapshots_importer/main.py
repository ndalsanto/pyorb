#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 10:24:47 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

#%%

import numpy as np

import sys
sys.path.insert(0, '../../core')

import rb_manager as rm
import affine_decomposition as ad

print( rm.__doc__ )






my_affine_decomposition = ad.AffineDecompositionHandler( )
my_affine_decomposition.set_Q( 4, 1 )                   # number of affine terms
my_affine_decomposition.import_affine_matrices( 'affine_matrix_A' )
my_affine_decomposition.import_affine_vectors(  'affine_vector_f' )


my_rb_manager = rm.RbManager( )
my_rb_manager.set_affine_decomposition_handler( my_affine_decomposition )

snapshots_file = 'train_snapshots_matrix_20_50.txt'
my_rb_manager.import_snapshots_matrix( snapshots_file )

my_rb_manager.build_rb_approximation( 10**(-4) )




my_rb_manager.print_rb_summary( )


