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

class fom_problem( ):

    def __init__( self, _parameter_handler ):

        self.define_theta_functions( )
        self.M_parameter_handler = _parameter_handler

        return

    def get_theta_a( self, _param, _q ):
        return self.M_theta_a( _param, _q )

    def get_theta_f( self, _param, _q ):
        return self.M_theta_f( _param, _q )

    def define_theta_functions( self ):

        em.error_raiser( 'SystemError', 'fom_problem::define_theta_functions', "You should define the theta function specific for your problem in the inherited class." )

        return

    # initialize anything which needs to be specified for using the external engine
    def configure_fom( self, _external_engine, _fom_specifics ):

        self.M_external_engine = _external_engine
        self.M_fom_specifics  = _fom_specifics
        self.M_configured_fom = True

        return

    def check_configured_fom( self ):

        if self.M_configured_fom == False:
            em.error_raiser( 'SystemError', 'fom_problem::retrieve_fom_data', "The fom problem has not been configured." )

        return

    def solve_fom_problem( self, _param ):
        self.check_configured_fom( )
        sol = self.M_external_engine.solve_parameter( _param, self.M_fom_specifics )
        return sol

    def compute_fom_product( self, _basis, _q, _operator ):

        print( "Performing compute_fom_product" )
        product = self.M_external_engine.build_rb_affine_component( _basis, _q, _operator, self.M_fom_specifics )

        return product.array

    def retrieve_fom_affine_components( self, _operator, _num_affine_components ):
        return self.M_external_engine.build_fom_affine_components( _operator, _num_affine_components, self.M_fom_specifics )

    def assemble_fom_matrix( self, _param, _elements=[], _indices=[] ):
        return self.M_external_engine.assemble_fom_matrix( _param, self.M_fom_specifics, _elements, _indices )

    def assemble_fom_rhs( self, _param, _elements=[], _indices=[] ):
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
        return self.M_external_engine.find_mdeim_elements_fom_specifics( self.M_fom_specifics, _indices_mat )

    def find_deim_elements_fom_specifics( self, _indices ):
        return self.M_external_engine.find_deim_elements_fom_specifics( self.M_fom_specifics, _indices )

    M_parameter_handler = None
    M_configured_fom = False

    # engine used to perform offline computation relying on an external engine
    M_external_engine = None
    M_fom_specifics = 0

    # theta functions
    M_theta_a = default_theta_function
    M_theta_f = default_theta_function
    M_current_parameter = np.zeros( 0 )
