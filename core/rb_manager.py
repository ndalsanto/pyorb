#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 10:22:53 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch

This module allows to handle RB methods in python starting from a set of FE arrays imported from another library.
"""



import numpy as np

import affine_decomposition as ad


class RbManager( ):
    
    def __init__( self ):
        return

    def import_snapshots_matrix( self, _input_file ):
        self.M_snapshots_matrix = np.loadtxt( _input_file )
        self.M_ns = self.M_snapshots_matrix.shape[1]
        return 

    def import_test_snapshots_matrix( self, _input_file ):
        self.M_test_snapshots_matrix = np.loadtxt( _input_file )
        self.M_ns_test = self.M_test_snapshots_matrix.shape[1]
        return 
    
    def get_snapshots_matrix( self ):
        return self.M_snapshots_matrix
    
    def get_number_of_snapshots( self ):
        return self.M_ns
    
    def get_basis( self ):
        return self.M_basis
    
    def get_number_of_basis( self ):
        return self.M_N
    
    def print_snapshots_matrix( self ):
        
        print( "The snapshots matrix is: \n" )
        print( self.M_snapshots_matrix )
        
        return
    
    def set_random_snapshots_matrix( self, _rows = 4, _cols = 2 ):

        self.M_snapshots_matrix = np.random.randn( _rows, _cols )
        
        return 
    
    def perform_pod( self, _tol = 10**(-5) ):
        
        U, s, V = np.linalg.svd( self.M_snapshots_matrix, full_matrices=False )
        
        total_energy = np.dot( s, np.transpose(s) )
        
        print( "The total energy of the field is %g" % total_energy )
        
        self.M_N = 0
        
        cumulative_energy = 0.0
        record_cumulative_energy = np.ones( self.M_ns )
        
        while cumulative_energy / total_energy < 1 - _tol**2:
            print( "Now N is %d and the cumulative energy is %g --- %g --- current sv %g " % ( self.M_N, cumulative_energy, cumulative_energy / total_energy, s[ self.M_N ] ) )
            record_cumulative_energy[ self.M_N ] = cumulative_energy
            cumulative_energy = cumulative_energy + s[ self.M_N ]*s[ self.M_N ]  # add the energy of next basis
            self.M_N = self.M_N + 1                                       # add a basis in the count
        
        print( "The number of selected RB functions is %d" % self.M_N )

        self.M_basis = U[:, 0:self.M_N]
            
        return

    def print_rb_summary( self ):
        
        print( "\n\n ------------  RB SUMMARY  ------------\n\n" )
        print( "Number of snapshots                    %d" % self.M_ns )
        print( "Number of selected RB functions        %d" % self.M_N )
        
        self.M_affineDecomposition.print_ad_summary( )
        
        return

    def set_affine_decomposition_handler( self, _affineDecomposition ):
        self.M_affineDecomposition = _affineDecomposition
        return

    def build_rb_approximation( self, _tol = 10**(-5) ):
        self.perform_pod( _tol )
        
        self.M_affineDecomposition.build_rb_affine_decompositions( )
        
        return

    M_ns = 0
    M_snapshots_matrix = np.zeros( 0 )
    M_ns_test = 0
    M_test_snapshots_matrix = np.zeros( 0 )

    M_N = 0
    M_basis = np.zeros( 0 )

    M_affineDecomposition = ad.AffineDecompositionHandler( )




























