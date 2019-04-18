#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 12:02:21 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import pyorb_core.error_manager as em
import numpy as np

def default_theta_function( _param, _q ):

    em.error_raiser( 'SystemError', 'default_theta_function', "You are using the default theta function, please provide specific ones for your problem " )

    pass

def default_full_theta_function( _param ):

    em.error_raiser( 'SystemError', 'default_full_theta_function', "You are using the default full theta function, please provide specific ones for your problem " )

    pass

class fom_problem( ):

    def __init__( self, _parameter_handler, _external_engine = None, _fom_specifics = None ):

        if _external_engine is not None and _fom_specifics is not None:
            self.configure_fom( _external_engine, _fom_specifics )
        
        self.define_theta_functions( )
        self.M_parameter_handler = _parameter_handler

        return

    def get_theta_a( self, _param, _q ):
        return self.M_theta_a( _param, _q )

    def get_theta_f( self, _param, _q ):
        return self.M_theta_f( _param, _q )

    def get_full_theta_a( self, _param ):
        return self.M_full_theta_a( _param )

    def get_full_theta_f( self, _param ):
        return self.M_full_theta_f( _param )

    def define_theta_functions( self ):
        em.error_raiser( 'SystemError', 'fom_problem::define_theta_functions', "You should define the theta function specific for your problem in the inherited class." )
        return

    # initialize anything which needs to be specified for using the external engine
    def configure_fom( self, _external_engine, _fom_specifics ):

        self.M_external_engine = _external_engine
        self.set_fom_specifics( _fom_specifics )
        self.M_external_engine.initialize_fom_simulation( _fom_specifics )
        self.assemble_fom_natural_norm_matrix( self.M_fom_specifics )
        self.M_configured_fom = True

        return

    def assemble_fom_natural_norm_matrix( self, _fom_specifics ):
        
        self.check_configured_fom( )
        self.M_natural_norm_matrix = self.M_external_engine.assemble_fom_natural_norm_matrix( self.M_fom_specifics )

    def set_fom_specifics( self, _fom_specifics ):

        self.M_fom_specifics = _fom_specifics

        return

    def update_fom_specifics( self, _fom_specifics_update ):
        
#        self.M_fom_specifics.update( _fom_specifics_update )
        print( "Updating the fom specifics dictionary" )

        for key in _fom_specifics_update:
            self.M_fom_specifics[key] = _fom_specifics_update[key]
        
        return 

    def clear_fom_specifics( self, _fom_specifics_update ):
        
        print( "Clearing the fom specifics dictionary" )

        for key in _fom_specifics_update:
            self.M_fom_specifics.pop( key )
        
        return 

    def check_configured_fom( self ):

        if self.M_configured_fom == False:
            em.error_raiser( 'SystemError', 'fom_problem::retrieve_fom_data', "The fom problem has not been configured." )

        return

    def compute_natural_norm( self, _solution ):
        self.check_configured_fom( )
        sol = self.M_external_engine.compute_natural_norm( _solution, self.M_fom_specifics )
        
        return sol

    def solve_fom_problem( self, _param ):
        self.check_configured_fom( )
        sol = self.M_external_engine.solve_parameter( _param, self.M_fom_specifics )
        return sol

    def compute_fom_product( self, _basis, _q, _operator ):

        print( "Performing compute_fom_product" )
        product = self.M_external_engine.build_rb_affine_component( _basis, _q, _operator, self.M_fom_specifics )

        return product.array

    def retrieve_fom_affine_components( self, _operator, _num_affine_components ):
        self.check_configured_fom( )
        return self.M_external_engine.build_fom_affine_components( _operator, _num_affine_components, self.M_fom_specifics )

    def assemble_fom_matrix( self, _param, _elements=[], _indices=[] ):
        self.check_configured_fom( )

        return self.M_external_engine.assemble_fom_matrix( _param, self.M_fom_specifics, _elements, _indices )

    def assemble_fom_rhs( self, _param, _elements=[], _indices=[] ):
        self.check_configured_fom( )
        return self.M_external_engine.assemble_fom_rhs( _param, self.M_fom_specifics, _elements, _indices )

    def get_num_parameters( self ):
        return self.M_parameter_handler.get_num_parameters( )

    def generate_parameter( self ):
        return self.M_parameter_handler.generate_parameter( )

    def get_parameter( self ):
        self.M_current_parameter = self.M_parameter_handler.get_parameter( )
        return self.M_current_parameter

    def get_parameter_handler( self ):
        return self.M_parameter_handler

    def find_mdeim_elements_fom_specifics( self, _indices_mat ):
        self.check_configured_fom( )
        return self.M_external_engine.find_mdeim_elements_fom_specifics( self.M_fom_specifics, _indices_mat )

    def find_deim_elements_fom_specifics( self, _indices ):
        self.check_configured_fom( )
        return self.M_external_engine.find_deim_elements_fom_specifics( self.M_fom_specifics, _indices )

    M_parameter_handler = None
    M_configured_fom = False

    # engine used to perform offline computation relying on an external engine
    M_external_engine = None
    M_fom_specifics = None

    M_natural_norm_matrix = None

    # theta functions
    M_theta_a = default_theta_function
    M_theta_f = default_theta_function
    M_full_theta_a = default_full_theta_function
    M_full_theta_f = default_full_theta_function
    M_current_parameter = np.zeros( 0 )
    





