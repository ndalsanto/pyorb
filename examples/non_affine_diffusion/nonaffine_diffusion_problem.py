#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 18:43:23 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import pyorb_core.pde_problem.fom_problem as fp

def ndp_theta_f( _param, _q = 0 ):
    return 1.0

class ndp_theta_fs( ):
    def __init__( self ):
        return
    
    def ndp_theta_f( self, _param, _q ):
        return self.M_deim.compute_deim_theta_coefficients_q( _param, _q )
    
    def ndp_full_theta_f( self, _param ):
        return self.M_deim.compute_deim_theta_coefficients( _param )

    def set_deim( self, _deim ):
        self.M_deim = _deim

    M_deim = None


class ndp_theta_As( ):
    
    def __init__( self ):
        return
    
    def ndp_theta_a( self, _param, _q ):
        return self.M_mdeim.compute_theta_coefficients_q( _param, _q )

    def ndp_full_theta_a( self, _param ):
        return self.M_mdeim.compute_theta_coefficients( _param )

    def set_mdeim( self, _mdeim ):
        self.M_mdeim = _mdeim

    M_mdeim = None



class nonaffine_diffusion_problem( fp.fom_problem ):

    def __init__( self, _parameter_handler, _external_engine = None, _fom_specifics = None ):
        fp.fom_problem.__init__( self, _parameter_handler, _external_engine, _fom_specifics )
        return
   
    def set_deim( self, _deim ):
        self.M_ndp_theta_fs.set_deim( _deim )
   
    def set_mdeim( self, _mdeim ):
        self.M_ndp_theta_As.set_mdeim( _mdeim )

    def define_theta_functions( self ):
        
        self.M_theta_a = self.M_ndp_theta_As.ndp_theta_a
        self.M_full_theta_a = self.M_ndp_theta_As.ndp_full_theta_a
        
        if self.M_fom_specifics['use_nonhomogeneous_dirichlet'] == 'Y':
            self.M_theta_f = self.M_ndp_theta_fs.ndp_theta_f
            self.M_full_theta_f = self.M_ndp_theta_fs.ndp_full_theta_f
        else:
            self.M_theta_f = ndp_theta_f
            self.M_full_theta_f = self.M_ndp_theta_fs.ndp_full_theta_f
        
        return
    
    def generate_parameter( self ):
        min_mu = self.get_parameter_handler( ).get_min_parameters( )
        max_mu = self.get_parameter_handler( ).get_max_parameters( )
        
        super(nonaffine_diffusion_problem, self).generate_parameter( )
        self.M_current_parameter = super(nonaffine_diffusion_problem, self).get_parameter( )

        self.M_current_parameter[2] = ( self.M_current_parameter[2] - min_mu[2] )**2 / (max_mu[2] - min_mu[2]) \
                                      + min_mu[2]
    
    M_ndp_theta_As = ndp_theta_As( )

    M_ndp_theta_fs = ndp_theta_fs( )
        






