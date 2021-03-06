#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 28 15:08:31 2019

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

#%%

import numpy as np

import sys
sys.path.insert(0, '../../')
print(sys.path)

import pyorb_core.tpl_managers.external_engine_manager as mee
import pyorb_core.error_manager as em

matlab_library_path = 'path/to/MATLAB/library'

if matlab_library_path == 'path/to/MATLAB/library':
    em.error_raiser( 'SystemError', 'Please specify the path for your matlab library' )

# playing around with engine manager 
my_matlab_engine_manager = mee.external_engine_manager( 'matlab', matlab_library_path )
my_matlab_engine_manager.start_engine( )
my_matlab_external_engine = my_matlab_engine_manager.get_external_engine( )

import pyorb_core.pde_problem.parameter_handler as ph

mu0_min = 0.01; mu0_max = 1.
mu1_min = -1.;  mu1_max = 1.
mu2_min = -1.;  mu2_max = 1.

param_min = np.array([mu0_min, mu1_min, mu2_min])
param_max = np.array([mu0_max, mu1_max, mu2_max])
num_parameters = param_min.shape[0]

# preparing the parameter handler
my_parameter_handler = ph.Parameter_handler( )
my_parameter_handler.assign_parameters_bounds( param_min, param_max )

# define the fem problem
import elliptic_example_problem as ell_ex

my_ell_ex = ell_ex.elliptic_example_problem( my_parameter_handler )

fem_size = 20

fom_specifics = { 
        'number_of_elements': fem_size, 
        'polynomial_degree' : 'P1',
        'model'             : 'elliptic_example',
        'use_nonhomogeneous_dirichlet' : 'N',
        'mesh_name'         : str(fem_size) + 'x' + str(fem_size) }

my_ell_ex.configure_fom( my_matlab_external_engine, fom_specifics )


#%%

import pyorb_core.rb_library.affine_decomposition as ad

# defining the affine decomposition structure
my_affine_decomposition = ad.AffineDecompositionHandler( )
my_affine_decomposition.set_Q( 2, 2 )               # number of affine terms

# building the RB manager
import pyorb_core.rb_library.rb_manager as rm
print( rm.__doc__ )
my_rb_manager = rm.RbManager( my_affine_decomposition, my_ell_ex )

my_rb_manager.build_rb_approximation( 50, 10**(-5) )

#%%
# printing summary
my_rb_manager.print_rb_offline_summary( )
avg_error = my_rb_manager.test_rb_solver( 10 )

#%%

my_matlab_engine_manager.quit_engine( )
