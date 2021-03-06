#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 7 21:25:00 2018

@author: Luca Pegolotti
@email : luca.pegolotti@epfl.ch
"""

import numpy as np
from mpi4py import MPI

import sys
sys.path.insert(0, '../../')
print(sys.path)

import pyorb_core.tpl_managers.external_engine_manager as mee

cpp_library_folder = '/usr/scratch/dalsanto/EPFL/DeepLearning/LifeV/lifev-env/pyorb-lifev-api-build/libpyorb-lifev-api.so'

# playing around with engine manager
my_cpp_engine_manager = mee.external_engine_manager( 'cpp', cpp_library_folder )
my_cpp_engine_manager.start_engine( )
my_cpp_external_engine = my_cpp_engine_manager.get_external_engine( )

import pyorb_core.pde_problem.parameter_handler as ph

mu0_min = 0.5; mu0_max = 10.
mu1_min = 0; mu1_max = 30.

param_min = np.array([mu0_min, mu1_min])
param_max = np.array([mu0_max, mu1_max])
num_parameters = param_min.shape[0]

# preparing the parameter handler
my_parameter_handler = ph.Parameter_handler( )
my_parameter_handler.assign_parameters_bounds( param_min, param_max )

# define the fem problem
import affine_advection as aap

my_aap = aap.affine_advection_problem( my_parameter_handler )

fom_specifics = {
        'model'             : 'affine_advection', \
        'datafile_path'     : './data'}

my_aap.configure_fom( my_cpp_external_engine, fom_specifics )

#%%

export_solution = 1

if export_solution:
    my_parameter_handler.generate_parameter( )
    param = my_parameter_handler.get_parameter( )
    # param = np.array([0.19891, 90])
    # param = np.array([0.7, 90])
    my_cpp_external_engine.solve_parameter_and_export_solution( param, fom_specifics )

#%%

import pyorb_core.rb_library.affine_decomposition as ad

# defining the affine decomposition structure
my_affine_decomposition = ad.AffineDecompositionHandler( )
my_affine_decomposition.set_Q( 3, 1 )               # number of affine terms

# building the RB manager
import pyorb_core.rb_library.rb_manager as rm
print( rm.__doc__ )
my_rb_manager = rm.RbManager( my_affine_decomposition, my_aap )

SAVE_OFFLINE = 1

name_str = "affine_advection"

if SAVE_OFFLINE == 1:
    my_rb_manager.save_offline_structures( "offline_" + name_str + "/snapshots_" + name_str + '.txt', \
                                           "offline_" + name_str + "/basis_" + name_str + '.txt', \
                                           "offline_" + name_str + "/rb_affine_components_" + name_str, \
                                           "offline_" + name_str + "/offline_parameters.data" )

my_rb_manager.build_rb_approximation( 150, 10**(-6) )

#%%
# printing summary
my_rb_manager.print_rb_offline_summary( )

avg_error = my_rb_manager.test_rb_solver( 10 )

#%%

my_cpp_engine_manager.quit_engine( )
