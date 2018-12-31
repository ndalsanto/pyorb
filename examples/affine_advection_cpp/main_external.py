#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 10:24:47 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch

An example where the RB method is constructed by importing from files the snapshots and the affine decomposition of FE matrices and vectors

"""

import sys
sys.path.insert(0, '../../')

#%%

import numpy as np


print(sys.path)

#import rb_manager as rm
#import affine_decomposition as ad
#import parameter_handler as ph

import pyorb_core.rb_library.rb_manager as rm
import pyorb_core.rb_library.affine_decomposition as ad
import pyorb_core.pde_problem.parameter_handler as ph

import affine_advection as aap

print( rm.__doc__ )

mu0_min = 0.5; mu0_max = 10.
mu1_min = 0; mu1_max = 30.

param_min = np.array([mu0_min, mu1_min])
param_max = np.array([mu0_max, mu1_max])
num_parameters = param_min.shape[0]

# preparing the parameter handler
my_parameter_handler = ph.Parameter_handler( )
my_parameter_handler.assign_parameters_bounds( param_min, param_max )

my_parameter_handler.generate_parameter( )
my_parameter_handler.print_parameters( )

# define the fem problem
my_aap = aap.affine_advection_problem( my_parameter_handler )

output_fdr = 'offline_affine_advection/'

# defining the affine decomposition structure
my_affine_decomposition = ad.AffineDecompositionHandler( )
my_affine_decomposition.set_Q( 3, 1 )                   # number of affine terms
my_affine_decomposition.import_rb_affine_matrices( output_fdr + 'rb_affine_components_affine_advection_A' )
my_affine_decomposition.import_rb_affine_vectors(  output_fdr + 'rb_affine_components_affine_advection_f' )

# building the RB manager
my_rb_manager = rm.RbManager( my_affine_decomposition, my_aap )

# importing snapshots, offline parameters and building RB space
my_rb_manager.import_snapshots_parameters( output_fdr + 'offline_parameters.data' )

snapshots_file = output_fdr + 'snapshots_affine_advection.txt'

my_rb_manager.import_snapshots_matrix( snapshots_file )

n_snapshots = 50
if my_rb_manager.M_ns < n_snapshots :
    print( 'We miss some snapshots! I have only %d in memory and I need to compute %d more.' % (my_rb_manager.M_ns, n_snapshots-my_rb_manager.M_ns) )
    my_rb_manager.build_snapshots( n_snapshots - my_rb_manager.M_ns )

my_rb_manager.perform_pod( 10**(-6) )

# printing summary
my_rb_manager.print_rb_offline_summary( )

my_rb_manager.import_test_parameters( output_fdr + 'offline_parameters.data' )
my_rb_manager.import_test_snapshots_matrix( snapshots_file )
for snapshot_number in range(100):
    my_rb_manager.compute_rb_test_snapshots_error( snapshot_number )
