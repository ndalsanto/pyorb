#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 31 18:43:23 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import pyorb_core.pde_problem.fom_problem as fp
import numpy as np

def ns_theta_f( _param, _q = 0 ):
    return 1.0

def ns_full_theta_f( _param ):
    return np.array([1.0])

class ns_theta_fs( ):
    def __init__( self ):
        return
    
    def ns_theta_f( self, _param, _q ):
        return self.M_deim.compute_deim_theta_coefficients_q( _param, _q )
    
    def ns_full_theta_f( self, _param ):
        return self.M_deim.compute_deim_theta_coefficients( _param )

    def set_deim( self, _deim ):
        self.M_deim = _deim

    M_deim = None


class ns_theta_As( ):
    
    def __init__( self ):
        return
    
    def ns_theta_a( self, _param, _q ):
        return self.M_mdeim.compute_theta_coefficients_q( _param, _q )

    def ns_full_theta_a( self, _param ):
        return self.M_mdeim.compute_theta_coefficients( _param )

    def set_mdeim( self, _mdeim ):
        self.M_mdeim = _mdeim

    M_mdeim = None



class navier_stokes_problem( fp.fom_problem ):

    def __init__( self, _parameter_handler, _external_engine = None, _fom_specifics = None ):
        fp.fom_problem.__init__( self, _parameter_handler, _external_engine, _fom_specifics )
        return
   
    def set_deim( self, _deim ):
        self.M_ns_theta_fs.set_deim( _deim )
   
    def set_mdeim( self, _mdeim ):
        self.M_ns_theta_As.set_mdeim( _mdeim )

    def define_theta_functions( self ):
        
#        self.M_theta_a = self.M_ns_theta_As.ns_theta_a
#        self.M_full_theta_a = self.M_ns_theta_As.ns_full_theta_a
#        
#        if self.M_fom_specifics['use_nonhomogeneous_dirichlet'] == 'Y':
#            self.M_theta_f = self.M_ns_theta_fs.ns_theta_f
#            self.M_full_theta_f = self.M_ns_theta_fs.ns_full_theta_f
#        else:
#            self.M_theta_f = ns_theta_f
#            self.M_full_theta_f = ns_full_theta_f
        
        return
    
    M_ns_theta_As = ns_theta_As( )
    M_ns_theta_fs = ns_theta_fs( )
        






