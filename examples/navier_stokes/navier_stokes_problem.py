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

        n_deim = self.M_deim.get_num_mdeim_basis( )

        if _q >= 0 and _q < n_deim:
            return self.M_deim.compute_deim_theta_coefficients_q( _param, _q )
        
        if _q == n_deim + 1:
            return 1.0

        # diffusion affine component
        if _q == n_deim:
            return 0.01 * _param[0]
    
    
    def ns_full_theta_f( self, _param ):
        return self.M_deim.compute_deim_theta_coefficients( _param )

    def set_deim( self, _deim ):
        self.M_deim = _deim

    M_deim = None


class ns_theta_As( ):
    
    def __init__( self ):
        return
    
    def ns_theta_a( self, _param, _q ):
        
        n_m_deim = self.M_mdeim.get_num_mdeim_basis( )
        
        if _q >= 0 and _q < n_m_deim:
            return self.M_mdeim.compute_theta_coefficients_q( _param, _q )

        if _q == n_m_deim:
            return 1.0

        # diffusion affine component
        if _q == n_m_deim + 1:
            return 0.01 * _param[0]

        

    def ns_full_theta_a( self, _param ):
        
        n_m_deim = self.M_mdeim.get_num_mdeim_basis( )
        thetas = np.zeros( n_m_deim + 2 )
        
        thetas[0:n_m_deim] = self.M_mdeim.compute_theta_coefficients( _param )
        thetas[n_m_deim] = _param[0]
        thetas[n_m_deim+1] = 1.0

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
        
        self.M_theta_a = self.M_ns_theta_As.ns_theta_a
        self.M_full_theta_a = self.M_ns_theta_As.ns_full_theta_a
        
        self.M_theta_f = ns_theta_f
        self.M_full_theta_f = ns_full_theta_f
        
        return
    
    def build_ns_rb_jacobian( self, _uN, _rb_affine_decomposition, _theta_a ):
        
        ns_jac = np.zeros( _rb_affine_decomposition.get_rb_affine_matrix(0).shape )

        for iQa in range( len(_theta_a) ):
            ns_jac = ns_jac + _theta_a[iQa] * _rb_affine_decomposition.get_rb_affine_matrix( iQa )

        for iQn in range( len(_uN) ):
            ns_jac = ns_jac + _uN[iQn] * _rb_affine_decomposition.get_rb_affine_matrix( iQn + len(_theta_a) )

        # the 2nd time is to include the flipped terms coming from the derivation
        for iQn in range( len(_uN) ):
            ns_jac = ns_jac + _uN[iQn] * _rb_affine_decomposition.get_rb_affine_matrix( iQn + len(_theta_a) + len(_uN) )
        
        return ns_jac

    def build_ns_rb_matrix( self, _uN, _rb_affine_decomposition, _theta_a ):
        
        ns_rb_mat = np.zeros( _rb_affine_decomposition.get_rb_affine_matrix(0).shape )

        for iQa in range( len(_theta_a) ):
            ns_rb_mat = ns_rb_mat + _theta_a[iQa] * _rb_affine_decomposition.get_rb_affine_matrix( iQa )
        
        for iQn in range( len(_uN) ):
            ns_rb_mat = ns_rb_mat + _uN[iQn] * _rb_affine_decomposition.get_rb_affine_matrix( iQn + len(_theta_a) )

        return ns_rb_mat

    def build_ns_rb_vector( self, _theta_f, _rb_affine_decomposition ):
        
        ns_rb_rhs = np.zeros( _rb_affine_decomposition.get_rb_affine_vector(0).shape )

        for iQf in range( len(_theta_f) ):
            ns_rb_rhs = ns_rb_rhs + _theta_f[iQf] * _rb_affine_decomposition.get_rb_affine_vector( iQf )

        return ns_rb_rhs

    def rb_residual( self, _uN, _rb_affine_decomposition, _theta_f, _theta_a ):

        ns_rb_mat = self.build_ns_rb_matrix( _uN, _rb_affine_decomposition, _theta_a )
        ns_rb_rhs = self.build_ns_rb_vector( _theta_f, _rb_affine_decomposition )

        res = ns_rb_mat.dot( _uN ) - ns_rb_rhs
        
        return res
    
    
    
    
    M_ns_theta_As = ns_theta_As( )
    M_ns_theta_fs = ns_theta_fs( )
        






