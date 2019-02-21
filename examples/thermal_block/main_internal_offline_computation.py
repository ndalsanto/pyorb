    #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:31:17 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch

An example where the RB method is constructed by solving the fem problem with MATLAB to compute the snapshots and the affine decomposition of FE matrices and vectors 

The affine decomposition is computed with MDEIM

"""

#%%

import numpy as np

import sys
sys.path.insert(0, '../../')
print(sys.path)

import pyorb_core.tpl_managers.external_engine_manager as mee

matlab_library_path = '/usr/scratch/dalsanto/EPFL/DeepLearning/feamat/'

# playing around with engine manager 
my_matlab_engine_manager = mee.external_engine_manager( 'matlab', matlab_library_path )
my_matlab_engine_manager.start_engine( )
my_matlab_external_engine = my_matlab_engine_manager.get_external_engine( )

import pyorb_core.pde_problem.parameter_handler as ph

mu0_min = 1.0; mu0_max = 50.
mu1_min = 1.0; mu1_max = 50.
mu2_min = 1.0; mu2_max = 50.

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
        'number_of_elements': 10, 
        'polynomial_degree' : 'P1',
        'model'             : 'thermal_block',
        'use_nonhomogeneous_dirichlet' : 'N',
        'mesh_name'                    : 'cube10x10' }

my_tbp.configure_fom( my_matlab_external_engine, fom_specifics )


#%%

import pyorb_core.rb_library.affine_decomposition as ad

# defining the affine decomposition structure
my_affine_decomposition = ad.AffineDecompositionHandler( )
my_affine_decomposition.set_Q( 4, 1 )               # number of affine terms

# building the RB manager
import pyorb_core.rb_library.rb_manager as rm
print( rm.__doc__ )
my_rb_manager = rm.RbManager( my_tbp, my_affine_decomposition )

my_rb_manager.build_rb_approximation( 50, 10**(-6) )

#%%
# printing summary
my_rb_manager.print_rb_offline_summary( )

avg_error = my_rb_manager.test_rb_solver( 10 )

#%%

my_matlab_engine_manager.quit_engine( )
