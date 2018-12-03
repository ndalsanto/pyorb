#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 18:43:23 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import pyorb_core.pde_problem.fom_problem as fp

def ndp_theta_f( _param, _q ):
    return 1.0

class ndp_theta_As( ):
    
    def __init__( self ):
        return
    
    def ndp_theta_a( self, _param, _q ):
        return self.M_mdeim.compute_theta_coefficients_q( _param, _q )

    def set_mdeim( self, _mdeim ):
        self.M_mdeim = _mdeim
        self.M_theta_a = self.M_mdeim

    M_mdeim = None



class nonaffine_diffusion_problem( fp.fom_problem ):

    def __init__( self, _parameter_handler ):
        fp.fom_problem.__init__( self, _parameter_handler )
        return
   
    def set_mdeim( self, _mdeim ):
        self.M_ndp_theta_As.set_mdeim( _mdeim )

    def define_theta_functions( self ):
        
        self.M_theta_a = self.M_ndp_theta_As.ndp_theta_a
        self.M_theta_f = ndp_theta_f
        
        return
    
    def generate_parameter( self ):
        min_mu = self.get_parameter_handler( ).get_min_parameters( )
        max_mu = self.get_parameter_handler( ).get_max_parameters( )
        
        super(nonaffine_diffusion_problem, self).generate_parameter( )
        self.M_current_parameter = super(nonaffine_diffusion_problem, self).get_parameter( )

        self.M_current_parameter[2] = ( self.M_current_parameter[2] - min_mu[2] )**2 / (max_mu[2] - min_mu[2]) \
                                      + min_mu[2]
    
    M_ndp_theta_As = ndp_theta_As( )
        
