#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:31:17 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch

An example where the RB method is constructed by solving the fem problem with MATLAB to compute the snapshots and the affine decomposition of FE matrices and vectors 

"""

#%%

import numpy as np

import sys
sys.path.insert(0, '../../')
print(sys.path)

import pyorb_core.tpl_managers.external_engine_manager as mee

matlab_library_path = '/usr/scratch/dalsanto/EPFL/DeepLearning/feamat/'
matlab_pyorb_interface = '/usr/scratch/dalsanto/EPFL/DeepLearning/pyorb-matlab-api/'

my_matlab_engine_manager = mee.external_engine_manager( 'matlab', matlab_library_path, matlab_pyorb_interface )
my_matlab_engine_manager.start_engine( )
my_matlab_external_engine = my_matlab_engine_manager.get_external_engine( )

import pyorb_core.pde_problem.parameter_handler as ph

mu0_min = 0.4; mu0_max = 0.6
mu1_min = 0.4; mu1_max = 0.6
mu2_min = 0.1; mu2_max = 0.5

param_min = np.array([mu0_min, mu1_min, mu2_min])
param_max = np.array([mu0_max, mu1_max, mu2_max])
num_parameters = param_min.shape[0]

# preparing the parameter handler
my_parameter_handler = ph.Parameter_handler( )
my_parameter_handler.assign_parameters_bounds( param_min, param_max )

# define the fem problem
import nonaffine_diffusion_problem as ndp

fem_size = 10
fem_size_str = str( fem_size )

fom_specifics = { 
        'number_of_elements': fem_size, \
        'polynomial_degree' : 'P1', \
        'model': 'nonaffine', \
        'simulation_name'   : 'nonaffine' + fem_size_str,\
        'use_nonhomogeneous_dirichlet' : 'Y', \
        'mesh_name' : fem_size_str + 'x' + fem_size_str, \
        'full_path'         : '/usr/scratch/dalsanto/EPFL/DeepLearning/pyorb_development/examples/non_affine_diffusion/' }

my_ndp = ndp.nonaffine_diffusion_problem( my_parameter_handler, my_matlab_external_engine, fom_specifics )

import pyorb_core.rb_library.m_deim as m_deim
my_mdeim = m_deim.Mdeim( my_ndp )
my_mdeim.perform_mdeim( 100, 10**(-6) )
my_ndp.set_mdeim( my_mdeim )

num_f_affine_components = 1

if fom_specifics['use_nonhomogeneous_dirichlet'] == 'Y':
    my_deim = m_deim.Deim( my_ndp )
    my_deim.perform_deim( 100, 10**(-6) )
    my_deim.print_reduced_indices( )
    my_ndp.set_deim( my_deim )
    num_f_affine_components = my_deim.get_num_basis( )
    
print( 'Number of affine basis for the rhs is %d ' % num_f_affine_components  )

import pyorb_core.rb_library.affine_decomposition as ad

# defining the affine decomposition structure
my_affine_decomposition = ad.AffineDecompositionHandler( )
my_affine_decomposition.set_Q( my_mdeim.get_num_mdeim_basis(), num_f_affine_components )               # number of affine terms

# we externally set the affine components for A, the ones for f are handled in the solver
my_affine_decomposition.set_affine_a( my_mdeim.get_basis_list( ) )

if fom_specifics['use_nonhomogeneous_dirichlet'] == 'Y':
    my_affine_decomposition.set_affine_f( my_deim.get_deim_basis_list( ) )

import pyorb_core.rb_library.rb_manager as rm
print( rm.__doc__ )
my_rb_manager = rm.RbManager( my_ndp, my_affine_decomposition )

SAVE_OFFLINE = 0

if SAVE_OFFLINE == 1:
    my_rb_manager.save_offline_structures( "offline_" + fem_size_str + "/test_snapshots_" + fem_size_str + '.txt', \
                                           "offline_" + fem_size_str + "/basis_" + fem_size_str + '.txt', \
                                           "offline_" + fem_size_str + "/rb_affine_components_" + fem_size_str, \
                                           'offline_' + fem_size_str + '/test_offline_parameters.data' )

my_rb_manager.build_rb_approximation( 100, 10**(-6) )

# printing summary
my_rb_manager.print_rb_offline_summary( )

my_rb_manager.test_rb_solver( 20 )

#%%

my_matlab_engine_manager.quit_engine( )
