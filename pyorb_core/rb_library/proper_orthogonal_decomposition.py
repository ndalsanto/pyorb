#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 17:21:31 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""
import numpy as np


class ProperOrthogonalDecompostion( ):
    
    def __init__( self ):
        return
        
    def perform_pod( self, _snapshots_matrix, _tol = 10**(-5) ):

        self.M_snapshots_matrix = _snapshots_matrix
        self.M_ns = self.M_snapshots_matrix.shape[1]
    
        U, s, V = np.linalg.svd( self.M_snapshots_matrix, full_matrices=False )
        
        total_energy = np.dot( s, np.transpose(s) )
        
        print( "The total energy of the field is %g" % total_energy )
        
        self.M_N = 0
        
        cumulative_energy = 0.0
        record_cumulative_energy = np.ones( self.M_ns )
        
        while cumulative_energy / total_energy < 1. - _tol**2 and self.M_N < self.M_ns:
            print( "Now N is %d and the cumulative energy is %g --- %g --- current sv %g " \
                                  % ( self.M_N, cumulative_energy, cumulative_energy / total_energy, s[ self.M_N ] ) )

            record_cumulative_energy[ self.M_N ] = cumulative_energy
            cumulative_energy = cumulative_energy + s[ self.M_N ] * s[ self.M_N ]   # add the energy of next basis
            self.M_N = self.M_N + 1                                                 # add a basis in the count

        print( "Final N is %d" % self.M_N )

        self.M_basis = U[:, 0:self.M_N]

    def get_basis( self ):
        return self.M_basis

    M_snapshots_matrix = np.zeros( ( 0, 0 ) )
    M_basis = np.zeros( ( 0, 0 ) )
    M_N = 0
    M_ns = 0
    