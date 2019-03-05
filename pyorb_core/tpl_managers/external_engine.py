#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:17:29 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""
import pyorb_core.error_manager as em

class external_engine( ):

    def __init__( self, _engine_type, _library_path ):
        
        self.M_engine_type  = _engine_type
        self.M_library_path = _library_path
        
        return

    def start_engine( self ):
        self.start_specific_engine( )
        return
    
    def quit_engine( self ):
        self.quit_specific_engine( )
        return

    def start_specific_engine( self ):
        
        em.error_raiser( 'SystemError', 'external_engine::start_engine', "You are using the default start_specific_engine, \
                          please provide specific ones for your specific engine " )
        return

    def quit_specific_engine( self ):
        
        em.error_raiser( 'SystemError', 'external_engine::quit_engine', "You are using the default quit_specific_engine, \
                          please provide specific ones for your specific engine " )
        return

    def convert_parameter( self, _param ):

        em.error_raiser( 'SystemError', 'external_engine::convert_parameter', "You are using the default convert_parameter, \
                          please provide specific ones for your specific engine " )
        return

    def solve_parameter( self, _param, _fom_specifics ):

        em.error_raiser( 'SystemError', 'external_engine::solve_parameter', "You are using the default solve_parameter, \
                          please provide specific ones for your specific engine " )
        return

    def build_fom_affine_components( self, _operator, _num_affine_components, _fom_specifics ):
        
        em.error_raiser( 'SystemError', 'external_engine::build_rb_affine_component', "You are using the default build_fom_affine_components, \
                          please provide specific ones for your specific engine " )
        return

    def assemble_fom_matrix( self, _param, _fom_specifics, _elements, _indices ):
        
        em.error_raiser( 'SystemError', 'external_engine::assemble_fom_matrix', "You are using the default assemble_fom_matrix, \
                          please provide specific ones for your specific engine " )
        return

    def find_mdeim_elements_fom_specifics( self, _fom_specifics, _indices_mat ):

        em.error_raiser( 'SystemError', 'external_engine::find_mdeim_elements_fom_specifics', \
                         "You are using the default find_mdeim_elements_fom_specifics, \
                          please provide specific ones for your specific engine " )
        return


    M_engine_type = ""
    M_library_path = ""
    M_engine = None


