#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 12:07:21 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

#%%

import numpy as np

import sys
sys.path.insert(0, '../../core')
sys.path.insert(0, '../thermal_block')
print(sys.path)


import manage_external_engine as mee
# playing around with engine manager 
my_cpp_engine_manager = mee.external_engine_manager( 'cpp', \
                                                     '/usr/scratch/dalsanto/EPFL/DeepLearning/LifeV/lifev-env/lifev-epfl-build/rb/liblifevreducedbasis.so' )

my_cpp_engine_manager.start_engine( )
my_cpp_engine = my_cpp_engine_manager.get_external_engine( )


import parameter_handler as ph

mu0_min = 1.0; mu0_max = 10.
mu1_min = 1.0; mu1_max = 10.
mu2_min = 1.0; mu2_max = 10.0

param_min = np.array([mu0_min, mu1_min, mu2_min])
param_max = np.array([mu0_max, mu1_max, mu2_max])
num_parameters = param_min.shape[0]

# preparing the parameter handler
my_parameter_handler = ph.Parameter_handler( )
my_parameter_handler.assign_parameters_bounds( param_min, param_max )

# define the fem problem 
import thermal_block_problem as tbp

my_tbp = tbp.thermal_block_problem( my_parameter_handler )

fom_specifics = { 
        'datafile_path'             : 'lifev_data/data',
        'model'                     : 'thermal_block'
        }

my_tbp.configure_fom( my_cpp_engine, fom_specifics )

import affine_decomposition as ad

# defining the affine decomposition structure
my_affine_decomposition = ad.AffineDecompositionHandler( )
my_affine_decomposition.set_Q( 4, 1 )               # number of affine terms

#my_affine_decomposition.import_affine_matrices( 'AAA_' )
#my_affine_decomposition.import_affine_vectors(  'fff_' )

# building the RB manager
import rb_manager as rm
print( rm.__doc__ )
my_rb_manager = rm.RbManager( my_affine_decomposition, my_tbp )

my_rb_manager.build_rb_approximation( 50, 10**(-6) )

my_rb_manager.print_affine_components( )



avg_error = my_rb_manager.test_rb_solver( 10 )




#%%
import numpy as np

import sys
sys.path.insert(0, '../../core')
sys.path.insert(0, '../thermal_block')
print(sys.path)

import manage_external_engine as mee
# playing around with engine manager 
my_cpp_engine_manager = mee.external_engine_manager( 'cpp', \
                                                     '/usr/scratch/dalsanto/EPFL/DeepLearning/LifeV/lifev-env/lifev-epfl-build/rb/liblifevreducedbasis.so' )

my_cpp_engine_manager.start_engine( )
my_cpp_engine = my_cpp_engine_manager.get_external_engine( )

fom_specifics = { 
        'datafile_path'             : 'lifev_data/data',
        'model'                     : 'thermal_block'
        }

AAA = my_cpp_engine.build_fom_matrix_affine_components( 1, fom_specifics )

AAA['A0'][1, 0:2]







