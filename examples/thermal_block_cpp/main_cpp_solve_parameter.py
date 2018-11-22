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
sys.path.insert(0, '../../')
sys.path.insert(0, '../thermal_block')
print(sys.path)


import pyorb_core.tpl_managers.external_engine_manager as mee
# playing around with engine manager 
my_cpp_engine_manager = mee.external_engine_manager( 'cpp', \
                                                     '/usr/scratch/dalsanto/EPFL/DeepLearning/LifeV/lifev-env/lifev-epfl-build/rb/liblifevreducedbasis.so' )

my_cpp_engine_manager.start_engine( )
my_cpp_engine_manager = my_cpp_engine_manager.get_external_engine( )


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
        'datafile_path'             : 'lifev_data/data',
        'model'                     : 'thermal_block'
        }

my_tbp.configure_fom( my_cpp_engine_manager, fom_specifics )

my_parameter_handler.generate_parameter( )
param = my_parameter_handler.get_parameter( )
u = my_tbp.solve_fom_problem( param )
















