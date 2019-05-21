#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 17:21:31 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""
import numpy as np
import pyorb_core.algebraic_utils as alg_ut


class ProperOrthogonalDecompostion( ):
    
    def __init__( self ):
        return
        
    def perform_pod( self, _snapshots_matrix, _tol = 10**(-5) ):

        self.M_snapshots_matrix = _snapshots_matrix
        self.M_ns = self.M_snapshots_matrix.shape[1]
        
        U, self.M_s, self.M_V = np.linalg.svd( self.M_snapshots_matrix, full_matrices=False )
        
        total_energy = np.dot( self.M_s, np.transpose(self.M_s) )
        
        print( "The total energy of the field is %g" % total_energy )
        
        self.M_N = 0
        
        cumulative_energy = 0.0
        record_cumulative_energy = np.ones( self.M_ns )
        
        while cumulative_energy / total_energy < 1. - _tol**2 and self.M_N < self.M_ns:
            print( "Now N is %d and the cumulative energy is %g --- %g --- current sv %g " \
                                  % ( self.M_N, cumulative_energy, cumulative_energy / total_energy, self.M_s[ self.M_N ] ) )

            record_cumulative_energy[ self.M_N ] = cumulative_energy
            cumulative_energy = cumulative_energy + self.M_s[ self.M_N ] * self.M_s[ self.M_N ]   # add the energy of next basis
            self.M_N = self.M_N + 1                                                 # add a basis in the count

        print( "Final N is %d" % self.M_N )

        self.M_basis = U[:, 0:self.M_N]
        self.compute_snapshots_coefficient( )


    def perform_pod_matrix( self, _snapshots_matrix, _tol, _norm_matrix ):

        print( "Performing POD with matrix inner product !" )

        self.M_snapshots_matrix = _snapshots_matrix
        self.M_ns = self.M_snapshots_matrix.shape[1]
        
        C_matrix = np.zeros( (self.M_ns, self.M_ns) )

        for iNs in range( self.M_ns ):
            Auh = alg_ut.sparse_matrix_vector_mul( _norm_matrix, self.M_snapshots_matrix[:, iNs] )
            for jNs in range( self.M_ns ):
                C_matrix[ iNs, jNs ] = np.dot( Auh, self.M_snapshots_matrix[:, jNs] )

        self.M_eig_vals, self.M_eig_vecs = np.linalg.eig( C_matrix )
        self.M_eig_vals = self.M_eig_vals.real
        self.M_eig_vecs = self.M_eig_vecs.real

        print( 'Eigenvalues' )
        print( self.M_eig_vals )

        total_energy = np.sum( self.M_eig_vals )
        
        print( "The total energy of the field is %g" % total_energy )
        
        self.M_N = 0
        
        cumulative_energy = 0.0
        record_cumulative_energy = np.ones( self.M_ns )
        
        while cumulative_energy / total_energy < 1. - _tol**2 and self.M_N < self.M_ns:
            print( "Now N is %d and the cumulative energy is %g --- %g --- current eig_val %g " \
                                  % ( self.M_N, cumulative_energy, cumulative_energy / total_energy, self.M_eig_vals[ self.M_N ] ) )

            record_cumulative_energy[ self.M_N ] = cumulative_energy
            cumulative_energy = cumulative_energy + self.M_eig_vals[ self.M_N ]   # add the energy of next basis
            self.M_N = self.M_N + 1                                               # add a basis in the count

        print( "Final N is %d" % self.M_N )

        self.M_basis = self.M_snapshots_matrix.dot( self.M_eig_vecs[:, 0:self.M_N] ) \
                                                                    / np.sqrt( self.M_eig_vals[0:self.M_N] )

#        self.M_basis = self.M_snapshots_matrix.dot( self.M_eig_vecs ) / np.sqrt( self.M_eig_vals )
#        self.M_basis = self.M_basis[0:self.M_N]


    def get_basis( self ):
        return self.M_basis

    def compute_snapshots_coefficient( self ):
        
        self.M_snapshots_coefficients = np.zeros( self.M_V.shape )

        for ii in range( len(self.M_s) ):
            self.M_snapshots_coefficients[ii, :] = self.M_s[ii] * self.M_V[ii, :]

        return

    def get_snapshots_coefficient( self ):
        return self.M_snapshots_coefficients

    M_snapshots_matrix = np.zeros( ( 0, 0 ) )
    M_basis = np.zeros( ( 0, 0 ) )
    M_N = 0
    M_ns = 0
    M_s = None
    M_V = None
    M_snapshots_coefficients= None
    
    
    
    
    
    
    
    