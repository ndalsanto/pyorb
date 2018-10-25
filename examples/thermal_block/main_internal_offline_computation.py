#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:31:17 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np

import sys
sys.path.insert(0, '../../core')

import manage_external_engine as mee

# playing around with engine manager 
my_matlab_engine_manager = mee.external_engine_manager( 'matlab', '/usr/scratch/dalsanto/EPFL/DeepLearning/feamat' )
my_matlab_engine_manager.start_engine( )
my_matlab_engine = my_matlab_engine_manager.get_engine( )

#%%

import rb_manager as rm
import affine_decomposition as ad
import parameter_handler as ph

import thermal_block_problem as tbp

print( rm.__doc__ )

mu0_min = 1.0; mu0_max = 50.
mu1_min = 1.0; mu1_max = 50.
mu2_min = 1.0; mu2_max = 50.

param_min = np.array([mu0_min, mu1_min, mu2_min])
param_max = np.array([mu0_max, mu1_max, mu2_max])
num_parameters = param_min.shape[0]

# preparing the parameter handler
my_parameter_handler = ph.Parameter_handler( )
my_parameter_handler.assign_parameters_bounds( param_min, param_max )

my_parameter_handler.generate_parameter( )
my_parameter_handler.print_parameters( )

# define the fem problem 
my_tbp = tbp.thermal_block_problem( )

# defining the affine decomposition structure
my_affine_decomposition = ad.AffineDecompositionHandler( )
my_affine_decomposition.set_Q( 4, 1 )                   # number of affine terms
my_affine_decomposition.import_affine_matrices( 'affine_matrix_20_A' )
my_affine_decomposition.import_affine_vectors(  'affine_vector_20_f' )

# building the RB manager
my_rb_manager = rm.RbManager( my_affine_decomposition, my_tbp )

# importing snapshots, offline parameters and building RB space
my_rb_manager.import_snapshots_parameters( 'train_parameters.data' )

snapshots_file = 'train_snapshots_matrix_20_50.txt'

my_rb_manager.import_snapshots_matrix( snapshots_file )

my_rb_manager.set_save_basis_functions( True, "basis.txt" )

my_rb_manager.build_rb_approximation( 10**(-6) )



# printing summary
my_rb_manager.print_rb_offline_summary( )


my_rb_manager.import_test_parameters( 'test_parameters.data' )
my_rb_manager.import_test_snapshots_matrix( 'test_snapshots_matrix_20_20.txt' )


for snapshot_number in range(20):
    my_rb_manager.compute_rb_test_snapshots_error( snapshot_number )
    







#%%

my_matlab_engine_manager.quit_engine( )
