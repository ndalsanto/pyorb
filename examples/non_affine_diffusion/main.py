#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:31:17 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np

import sys
sys.path.insert(0, '../../core')
print(sys.path)


import manage_external_engine as mee

# playing around with engine manager 
my_matlab_engine_manager = mee.external_engine_manager( 'matlab', '/usr/scratch/dalsanto/EPFL/DeepLearning/feamat' )
my_matlab_engine_manager.start_engine( )
my_matlab_external_engine = my_matlab_engine_manager.get_external_engine( )


import parameter_handler as ph

mu0_min = 0.4; mu0_max = 0.6
mu1_min = 0.4; mu1_max = 0.6
mu2_min = 0.25; mu2_max = 0.5

param_min = np.array([mu0_min, mu1_min, mu2_min])
param_max = np.array([mu0_max, mu1_max, mu2_max])
num_parameters = param_min.shape[0]

# preparing the parameter handler
my_parameter_handler = ph.Parameter_handler( )
my_parameter_handler.assign_parameters_bounds( param_min, param_max )

# define the fem problem 
import nonaffine_diffusion_problem as ndp

my_ndp = ndp.nonaffine_diffusion_problem( my_parameter_handler )

fem_size = 10
fem_size_str = str( fem_size )

fom_specifics = { 
        'number_of_elements': fem_size, 
        'polynomial_degree' : 'P1',
        'model': 'nonaffine' }

my_ndp.configure_fom( my_matlab_external_engine, fom_specifics )

import m_deim
my_mdeim = m_deim.Mdeim( my_ndp )

my_mdeim.perform_mdeim( 100, 10**(-5) )

#%%

import affine_decomposition as ad

# defining the affine decomposition structure
my_affine_decomposition = ad.AffineDecompositionHandler( )
# my_affine_decomposition.set_Q( my_mdeim.get_num_mdeim_basis( ), 1 )               # number of affine terms

my_affine_decomposition.set_Q( my_mdeim.get_num_mdeim_basis(), 1 )               # number of affine terms

my_affine_decomposition.import_affine_matrices( 'mdeim_basis_A' + fem_size_str + '_' )
my_affine_decomposition.import_affine_vectors(  'affine_vector_' + fem_size_str + '_f' )

# building the RB manager
import rb_manager as rm
print( rm.__doc__ )
my_rb_manager = rm.RbManager( my_affine_decomposition, my_ndp )

#%%

my_rb_manager.save_offline_structures( "offline_" + fem_size_str + "/snapshots_" + fem_size_str + '.txt', \
                                       "offline_" + fem_size_str + "/basis_" + fem_size_str + '.txt', \
                                       "offline_" + fem_size_str + "/rb_affine_components_" + fem_size_str, \
                                       'offline_' + fem_size_str + '/offline_parameters.data' )
my_rb_manager.build_rb_approximation( 500, 10**(-5) )


my_rb_manager.save_offline_structures( "offline_" + fem_size_str + "/test_snapshots_" + fem_size_str + '.txt', \
                                       "offline_" + fem_size_str + "/basis_" + fem_size_str + '.txt', \
                                       "offline_" + fem_size_str + "/rb_affine_components_" + fem_size_str, \
                                       'offline_' + fem_size_str + '/test_offline_parameters.data' )
my_rb_manager.build_snapshots( 500 )



#%%
# printing summary
my_rb_manager.print_rb_offline_summary( )

# my_affine_decomposition.print_affine_components( )




"""
my_rb_manager.import_test_parameters( 'test_parameters.data' )
my_rb_manager.import_test_snapshots_matrix( 'test_snapshots_matrix_20_20.txt' )


for snapshot_number in range(20):
    my_rb_manager.compute_rb_test_snapshots_error( snapshot_number )
    
avg_error = my_rb_manager.test_rb_solver( 100 )

"""

#%%

my_matlab_engine_manager.quit_engine( )
