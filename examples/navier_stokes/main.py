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

# playing around with engine manager
my_matlab_engine_manager = mee.external_engine_manager( 'matlab', matlab_library_path )
my_matlab_engine_manager.start_engine( )
my_matlab_external_engine = my_matlab_engine_manager.get_external_engine( )


import pyorb_core.pde_problem.parameter_handler as ph

mu0_min = 1.; mu0_max = 10.
mu1_min = 0.; mu1_max = 0.4

param_min = np.array([mu0_min, mu1_min])
param_max = np.array([mu0_max, mu1_max])
num_parameters = param_min.shape[0]

# preparing the parameter handler
my_parameter_handler = ph.Parameter_handler( )
my_parameter_handler.assign_parameters_bounds( param_min, param_max )

# define the fem problem
import navier_stokes_problem as ns

mesh = 'fine'

fom_specifics = {
        'model': 'navier_stokes',
        'mesh_name' : '/usr/scratch/dalsanto/EPFL/DeepLearning/DLPDEs/elliptic_example/navier_stokes_2d/bifurcation_' + mesh + '.msh' }

my_ns = ns.navier_stokes_problem( my_parameter_handler, my_matlab_external_engine, fom_specifics )

my_ns.generate_parameter( )
param = my_ns.get_parameter( )

do_offline = 0

base_foler = base_foler + '/'

#%%
import pyorb_core.rb_library.m_deim as m_deim
my_mdeim = m_deim.Mdeim( my_ns )

ns_m_deim = 500 

if do_offline == 1:
    my_mdeim.set_save_offline( True, base_foler + '/' )
#    my_mdeim.perform_mdeim( ns_mdeim, 10**(-6) )

    my_mdeim.build_mdeim_snapshots( ns_m_deim ) 
    my_mdeim.build_deim_basis( 10**(-6) ) 
else:
    my_mdeim.load_mdeim_basis( base_foler + '/' )

#my_ns.set_mdeim( my_mdeim )



#%%
my_deim = m_deim.Deim( my_ns )

if do_offline == 1:
    my_deim.set_save_offline( True, base_foler + '/' )
#    my_mdeim.perform_mdeim( ns_mdeim, 10**(-6) )

    my_deim.build_deim_snapshots( ns_m_deim ) 
    my_deim.build_deim_basis( 10**(-6) ) 
else:
    my_deim.load_deim_basis( base_foler + '/' )

num_f_affine_components = my_deim.get_num_basis( )

#my_ns.set_deim( my_deim )

print( 'Number of affine basis for the rhs is %d ' % num_f_affine_components  )
#%% 

import pyorb_core.rb_library.rb_manager as rm
print( rm.__doc__ )
my_rb_manager = rm.RbManager( my_ns )

SAVE_OFFLINE = 1

if SAVE_OFFLINE == 1:
    my_rb_manager.save_offline_structures( base_foler + "snapshots_" + mesh + '.txt', \
                                           base_foler + "basis_" + mesh + '.txt', \
                                           base_foler + "rb_affine_components_" + mesh, \
                                           base_foler + 'offline_parameters.data' )

n_s = 200

if do_offline == 1:
    my_rb_manager.build_snapshots( n_s )
else:
    my_rb_manager.import_snapshots_matrix( base_foler + "snapshots_" + mesh + '.txt' )
    my_rb_manager.import_snapshots_parameters( base_foler + 'offline_parameters.data' )

my_rb_manager.perform_pod( 10**(-4) )

rb_functions_dict = my_rb_manager.get_rb_functions_dict( )
my_rb_manager.update_fom_specifics( rb_functions_dict )

import pyorb_core.rb_library.affine_decomposition as ad

# defining the affine decomposition structure
my_affine_decomposition = ad.AffineDecompositionHandler( )
my_affine_decomposition.set_Q( my_mdeim.get_num_mdeim_basis()+my_rb_manager.get_number_of_basis()+1, \
                               num_f_affine_components )               # number of affine terms

# we externally set the affine components for A, the ones for f are handled in the solver
my_affine_decomposition.set_affine_a( my_mdeim.get_basis_list( ) )
my_affine_decomposition.set_affine_f( my_deim.get_deim_basis_list( ) )

my_rb_manager.set_affine_decomposition_handler( my_affine_decomposition )

#%%

if do_offline == 1:
    my_rb_manager.build_rb_affine_decompositions( )
    
    if SAVE_OFFLINE == 1:
        my_rb_manager.save_rb_affine_decomposition( )
else:
    my_affine_decomposition.import_rb_affine_matrices( base_foler + 'rb_affine_components_' + mesh + '_A' )
    my_affine_decomposition.import_rb_affine_vectors(  base_foler + 'rb_affine_components_' + mesh + '_f' )






#%%

#my_matlab_engine_manager.quit_engine( )
