#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 31 14:31:17 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch

An example where the RB method is constructed by solving the fem problem with MATLAB to compute the snapshots and the affine decomposition of FE matrices and vectors 

"""

#%%


def compute_rb_errors( n_tests, my_rb_manager, my_ns ):

  errors = np.zeros( n_tests ) 
 
  for ii in range( n_tests ):
    print('Parameter %d ' % ii)
    my_ns.generate_parameter( )
    param_1 = my_ns.get_parameter( )
    print( param_1 )

    start = time.time()
    un = my_ns.solve_rb_ns_problem( param_1, my_rb_manager.M_affineDecomposition )
    end = time.time()
    time_to_solve = end - start
    print( 'Time to solve RB problem %f ' % time_to_solve )

    s1 = my_ns.solve_fom_problem( param_1 )
    utildeh = my_rb_manager.reconstruct_fem_solution( un )
    e_h = s1 - utildeh
    norm_error = np.linalg.norm(e_h) / np.linalg.norm(s1)
    print( 'Norm of error %f' % norm_error )

    err = my_ns.compute_natural_norm( e_h ) / my_ns.compute_natural_norm( s1 )
    print( 'Norm of H1 error %f' % err )

    errors[ii] = err

  print('Mean error %f ' % np.mean(errors))



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

mu0_min = 1.; mu0_max = 10.
mu1_min = 0.; mu1_max = 0.1

param_min = np.array([mu0_min, mu1_min])
param_max = np.array([mu0_max, mu1_max])
num_parameters = param_min.shape[0]

# preparing the parameter handler
my_parameter_handler = ph.Parameter_handler( )
my_parameter_handler.assign_parameters_bounds( param_min, param_max )
param_range = 'param_03_'

do_offline = 0

offline_selection = ''

if mu1_max <= 0.250001:

    param_range = 'param_025_'
    ns_m_deim = 500
    n_s = 750
    ns_test = 500
    rb_tol = 10**(-4)
    do_offline = 0

    # fine
    param_range = 'param_025_'
    ns_m_deim = 750
    n_s = 1500
    ns_test = 50
    rb_tol = 10**(-4)
    do_offline = 0
    mu0_grid = 30
    mu1_grid = 50
    offline_selection = 'tensor'

if mu1_max <= 0.150001:

    param_range = 'param_015_'
    ns_m_deim = 200
    n_s = 500
    ns_test = 10
    rb_tol = 10**(-5)
    do_offline = 0
    
#    #fine
#    param_range = 'param_015_'
#    ns_m_deim = 200
#    n_s = 750
#    ns_test = 100
#    rb_tol = 2.*10**(-4)
#    do_offline = 0

if mu1_max <= 0.101:

    param_range = 'param_010_'
    ns_m_deim = 30
    n_s = 42
    mu0_grid = 6
    mu1_grid = 7
    offline_selection = 'tensor'
    ns_test = 10
    rb_tol = 10**(-4)
    do_offline = 1

if mu1_max < 0.0999:
    param_range = 'small_param_'
    ns_m_deim = 10 
    n_s = 40
    ns_test = 10
    rb_tol = 10**(-6)
    do_offline = 0

if mu1_max < 0.001:
    param_range = 'affine_param_'
    ns_m_deim = 3
    n_s = 300
    ns_test = 40
    rb_tol = 10**(-6)
    do_offline = 0

# define the fem problem
import navier_stokes_problem as ns

mesh = 'very_coarse'
#mesh = 'fine'

fom_specifics = {
        'model': 'navier_stokes', \
        'simulation_name'   : 'navier_stokes_' + mesh,\
        'use_nonhomogeneous_dirichlet' : 'Y', \
        'mesh_name' : '/usr/scratch/dalsanto/EPFL/DeepLearning/pyorb_development/examples/navier_stokes/meshes/bifurcation_' + mesh + '.msh', \
        'full_path'         : '/usr/scratch/dalsanto/EPFL/DeepLearning/pyorb_development/examples/navier_stokes/' }

my_ns = ns.navier_stokes_problem( my_parameter_handler, my_matlab_external_engine, fom_specifics )

##%%

my_ns.generate_parameter( )
param = my_ns.get_parameter( )

base_folder = 'offline_' + param_range + offline_selection + '_' + mesh + '/'
print(base_folder)

#%%
import pyorb_core.utils.array_utils as pyorb_array_utils
import pyorb_core.rb_library.m_deim as m_deim
my_mdeim = m_deim.Mdeim( my_ns )
my_deim  = m_deim.Deim(  my_ns )

if do_offline == 1:
    my_mdeim.set_save_offline( True, base_folder + '/' )
    my_mdeim.perform_mdeim( ns_m_deim, 10**(-6) )

    theta_mdeim_min, theta_mdeim_max = my_mdeim.compute_theta_min_max( )
    pyorb_array_utils.save_vector( theta_mdeim_min, base_folder + 'theta_mdeim_min.txt' )
    pyorb_array_utils.save_vector( theta_mdeim_max, base_folder + 'theta_mdeim_max.txt' )

    pyorb_array_utils.save_matrix( my_mdeim.M_snapshots_matrix, base_folder + 'matrix_snapshots.txt' )
    pyorb_array_utils.save_matrix( my_mdeim.M_snapshots_coefficients, base_folder + 'matrix_snapshots_theta.txt' )

else:
    my_mdeim.load_mdeim_offline( base_folder )
    theta_mdeim_min = np.loadtxt( base_folder + 'theta_mdeim_min.txt' )
    theta_mdeim_max = np.loadtxt( base_folder + 'theta_mdeim_max.txt' )
    my_mdeim.M_snapshots_coefficients = np.loadtxt( base_folder + 'matrix_snapshots_theta.txt' )
    my_mdeim.M_snapshots_matrix = np.loadtxt( base_folder + 'matrix_snapshots.txt' )

my_mdeim.M_snapshots_coefficients
num_mdeim_affine_components_A = my_mdeim.get_num_mdeim_basis( )

if do_offline == 1:
    my_deim.set_save_offline( True, base_folder + '/' )
    my_deim.perform_deim( ns_m_deim, 10**(-6) )
    theta_deim_min, theta_deim_max = my_deim.compute_theta_min_max( )
    pyorb_array_utils.save_vector( theta_deim_min, base_folder + 'theta_deim_min.txt' )
    pyorb_array_utils.save_vector( theta_deim_max, base_folder + 'theta_deim_max.txt' )

    pyorb_array_utils.save_matrix( my_deim.M_snapshots_matrix, base_folder + 'vector_snapshots.txt' )
    pyorb_array_utils.save_matrix( my_deim.M_snapshots_coefficients, base_folder + 'vector_snapshots_theta.txt' )

else:
    my_deim.load_deim_offline( base_folder + '/' )
    theta_deim_min = np.loadtxt( base_folder + 'theta_deim_min.txt' )
    theta_deim_max = np.loadtxt( base_folder + 'theta_deim_max.txt' )
    my_deim.M_snapshots_coefficients = np.loadtxt( base_folder + 'vector_snapshots_theta.txt' )
    my_deim.M_snapshots_matrix = np.loadtxt( base_folder + 'vector_snapshots.txt' )
    
my_deim.M_snapshots_coefficients
num_deim_affine_components_f = my_deim.get_num_basis( )
    
my_ns.set_mdeim( my_mdeim )
my_ns.set_deim( my_deim )

num_f_affine_components = my_deim.get_num_basis( )

#my_ns.set_deim( my_deim )

print( 'Number of affine basis for the rhs is %d ' % num_f_affine_components  )
##%%

import pyorb_core.rb_library.rb_manager as rm
print( rm.__doc__ )
my_rb_manager = rm.RbManager( my_ns )

SAVE_OFFLINE = 1

if SAVE_OFFLINE == 1:
    my_rb_manager.save_offline_structures( base_folder + "snapshots_" + mesh + '.txt', \
                                           base_folder + "basis_" + mesh + '.txt', \
                                           base_folder + "rb_affine_components_" + mesh, \
                                           base_folder + 'offline_parameters.data' )

if do_offline == 1:
    my_rb_manager.build_snapshots( n_s )
else:
    my_rb_manager.import_snapshots_matrix( base_folder + "snapshots_" + mesh + '.txt' )
    my_rb_manager.import_snapshots_parameters( base_folder + 'offline_parameters.data' )

my_rb_manager.perform_pod( rb_tol )

import pyorb_core.rb_library.affine_decomposition as ad

# want_to_solve_newton should be 1 to consider the affine components for the Jacobian
want_to_solve_newton = 0
want_to_solve_newton = 1

# defining the affine decomposition structure
# the two + 1 in the affine decoimpositions are the diffusion term and the lifting contribution in the nonlinear term
my_affine_decomposition = ad.AffineDecompositionHandler( )
my_affine_decomposition.set_Q( num_mdeim_affine_components_A + (1+want_to_solve_newton)*my_rb_manager.get_number_of_basis()+1+1, \
                               num_f_affine_components+1+1 )               # number of affine terms

# we externally set the affine components for A, the ones for f are handled in the solver
my_affine_decomposition.set_affine_a( my_mdeim.get_basis_list( ) )
my_affine_decomposition.set_affine_f( my_deim.get_deim_basis_list( ) )

my_rb_manager.set_affine_decomposition_handler( my_affine_decomposition )

#do_offline = 0

if do_offline == 1:
    rb_functions_dict = my_rb_manager.get_rb_functions_dict( )
    my_ns.update_fom_specifics( rb_functions_dict )

    my_rb_manager.build_rb_affine_decompositions( _build_rb_tpl=True )
    
    if SAVE_OFFLINE == 1:
        my_rb_manager.save_rb_affine_decomposition( )
    
#    clearing the fom specifics dictionary
    my_ns.clear_fom_specifics( rb_functions_dict )

else:
    my_affine_decomposition.import_rb_affine_matrices( base_folder + 'rb_affine_components_' + mesh + '_A' )
    my_affine_decomposition.import_rb_affine_vectors(  base_folder + 'rb_affine_components_' + mesh + '_f' )

my_rb_manager.print_rb_offline_summary( )

#%%

print( 'Solving RB problem ' )

train_parameters = my_rb_manager.get_offline_parameters( )

param_1 = train_parameters[0, :]
import time
start = time.time()
un = my_ns.solve_rb_ns_problem( param_1, my_rb_manager.M_affineDecomposition )
end = time.time()
time_to_solve = end - start
print( 'Time to solve RB problem %f ' % time_to_solve )

print('un is')
print(un)

s1 = my_rb_manager.M_snapshots_matrix[:, 0]
utildeh = my_rb_manager.reconstruct_fem_solution( un )
e_h = s1 - utildeh
norm_error = np.linalg.norm(e_h) / np.linalg.norm(s1)
print( 'Norm of error %f' % norm_error )

err = my_ns.compute_natural_norm( e_h ) / my_ns.compute_natural_norm( s1 )
print( err )

compute_rb_errors( 20, my_rb_manager, my_ns )


my_matlab_engine_manager.quit_engine( )
