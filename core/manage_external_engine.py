#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:17:29 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import matlab.engine
import error_manager as em

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

    def build_rb_affine_component( self, _param, _fom_specifics ):

        em.error_raiser( 'SystemError', 'external_engine::build_rb_affine_component', "You are using the default build_rb_affine_component, \
                          please provide specific ones for your specific engine " )
        return

    def build_fe_affine_components( self, _operator, _fom_specifics ):
        
        em.error_raiser( 'SystemError', 'external_engine::build_rb_affine_component', "You are using the default build_fe_affine_components, \
                          please provide specific ones for your specific engine " )
        return

    def assemble_fom_matrix( self, _param, _fom_specifics ):
        
        em.error_raiser( 'SystemError', 'external_engine::assemble_fom_matrix', "You are using the default assemble_fom_matrix, \
                          please provide specific ones for your specific engine " )
        return

    def find_mdeim_elements_fem_specifics( self, _fom_specifics, _indices_mat ):

        em.error_raiser( 'SystemError', 'external_engine::find_mdeim_elements_fem_specifics', \
                         "You are using the default find_mdeim_elements_fem_specifics, \
                          please provide specific ones for your specific engine " )
        return

        
    M_engine_type = ""
    M_library_path = ""
    M_engine = 0


class matlab_external_engine( external_engine ):
    
    def __init__( self, _engine_type, _library_path ):
        
        external_engine.__init__( self, _engine_type, _library_path )
        
        return

    def start_specific_engine( self ):
        
        self.M_engine = matlab.engine.start_matlab( )
        self.M_engine.addpath( self.M_engine.genpath( self.M_library_path ) )

        print( 'Successfully started matlab engine and corresponding FOM library %s ' % self.M_library_path )

        return

    def quit_specific_engine( self ):

        self.M_engine.quit( )

        print( 'Successfully quitted matlab engine' )

        return

    def convert_parameter( self, _param ):

        return matlab.double(_param.tolist())

    def convert_indices( self, _indices ):

        return matlab.int64(_indices.tolist())

    def solve_parameter( self, _param, _fom_specifics ):

        return self.M_engine.solve_parameter( self.convert_parameter( _param ), _fom_specifics )

    def build_rb_affine_component( self, _basis, _q, _operator, _fom_specifics ):

        return self.M_engine.build_rb_affine_component( _basis, _q, _operator, _fom_specifics )

    def build_fe_affine_components( self, _operator, _fom_specifics ):
        
        return self.M_engine.build_fe_affine_components( _operator, _fom_specifics )

    def assemble_fom_matrix( self, _param, _fom_specifics ):
        
        return self.M_engine.assemble_fom_matrix( self.convert_parameter( _param ), _fom_specifics )
    
    def find_mdeim_elements_fem_specifics( self, _fom_specifics, _indices_mat ):

        return self.M_engine.find_mdeim_elements_fem_specifics( _fom_specifics, self.convert_indices( _indices_mat ) )



class external_engine_manager( ):
    
    def __init__( self, _engine_type, _library_path ):
        
        self.M_engine_type  = _engine_type
        self.M_library_path = _library_path
        
        if _engine_type == 'matlab':
            self.M_external_engine = matlab_external_engine( _engine_type, _library_path )
        
        return

    M_engine_type = ""
    M_library_path = ""
    M_external_engine = 0
    
    def get_external_engine( self ):
        return self.M_external_engine

    def start_engine( self ):
        self.M_external_engine.start_engine( )
        return
    
    def quit_engine( self ):
        self.M_external_engine.quit_engine( )
        return



