#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 11 12:02:21 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

def default_theta_function( _param, _q ):

    print( "You are using the default theta function, please provide specific ones for your problem " )
    
    pass

class fom_problem( ):
    
    def __init__( self ):
    
        self.define_theta_functions( )
        
        return
    
    def get_theta_a( self, _param, _q ):
        return self.M_theta_a( _param, _q )
    
    def get_theta_f( self, _param, _q ):
        return self.M_theta_f( _param, _q )
    
    def define_theta_functions( self ):
        
        print( "You should define the theta function specific for your problem in the inherited class" )
        
        return
    
    # initialize anything which needs to be specified for using the external engine 
    def configure_fom( self, _externl_engine, _fom_specifics ):
        
        self.M_externl_engine = _externl_engine
        self.M_fom_specifics  = _fom_specifics
        self.M_configured_fom = True
        return

    # used to retrieve information from fom (e.g. number of dofs, fields, ... )
    def retrieve_fom_data( ):
        
        self.M_externl_engine.
        
        return
    
    def solve_fom_problem( self, _param ):
        
        self.solve_external_fom_problem( _param )
        
        return
    
    def solve_external_fom_problem( self, _param ):
        
        print( "You should define the solve_param function specific for your problem \
                and the external library you are using" )

        return

    def convert_parameter( self, _param ):
        return self.M_externl_engine( _param )


    M_configured_fom = False

    # engine used to perform offline computation relying on an external engine
    M_externl_engine = 0
    M_fom_specifics = 0
    
    # theta functions
    M_theta_a = default_theta_function
    M_theta_f = default_theta_function
    
    
    
    
    
    
    
    
    
    