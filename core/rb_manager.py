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
import fem_problem as fm


class RbHandler( ):
    
    def __init__( self ):
        return




class RbManager( ):
    
    def __init__( self, _affine_decomposition, _fem_problem ):
        
        self.set_affine_decomposition_handler( _affine_decomposition )
        self.set_fem_problem( _fem_problem )
        return

    def import_snapshots_parameters( self, _input_file ):
        self.M_offline_ns_parameters = np.loadtxt( _input_file )
        return 

    def import_snapshots_matrix( self, _input_file ):
        self.M_snapshots_matrix = np.loadtxt( _input_file )
        self.M_ns = self.M_snapshots_matrix.shape[1]
        return 

    def import_test_parameters( self, _input_file ):
        self.M_test_parameters = np.loadtxt( _input_file )
        return

    def import_test_snapshots_matrix( self, _input_file ):
        self.M_test_snapshots_matrix = np.loadtxt( _input_file )
        self.M_ns_test = self.M_test_snapshots_matrix.shape[1]
        return 

    def get_offline_parameter( self, _iP ):
        return self.M_offline_ns_parameters[_iP, :]
    
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
    
    def set_save_basis_functions( self, _save, _file ):
        self.M_save_basis_functions = _save
        self.M_save_file_basis_functions = _file
    
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

        self.M_basis = U[:, 0:self.M_N]
    
        if self.M_save_basis_functions == True:
            output_file = open( self.M_save_file_basis_functions, 'w+' )
                
            for iNs in range( self.M_basis.shape[0] ):
                for iP in range( self.M_basis.shape[1] ):
                    output_file.write( "%.10g" % self.M_basis[iNs, iP] )
    
                    if iP < self.M_basis.shape[1] - 1:
                        output_file.write( " " % self.M_basis[iNs, iP] )
                    else:
                        output_file.write( "\n" % self.M_basis[iNs, iP] )
        
            output_file.close( )

        return

    def print_rb_offline_summary( self ):
        
        print( "\n\n ------------  RB SUMMARY  ------------\n\n" )
        print( "Number of snapshots                    %d" % self.M_ns )
        print( "Number of selected RB functions        %d" % self.M_N )
        
        self.M_affineDecomposition.print_ad_summary( )
        
        return

    def set_affine_decomposition_handler( self, _affineDecomposition ):
        self.M_affineDecomposition = _affineDecomposition
        return

    def set_fem_problem( self, _fem_problem ):
        self.M_fem_problem = _fem_problem
        return


    def build_rb_approximation( self, _tol = 10**(-5) ):
        
        self.perform_pod( _tol )
        
        self.M_affineDecomposition.build_rb_affine_decompositions( self.M_basis )
        
        return

    M_ns = 0
    M_snapshots_matrix = np.zeros( 0 )
    M_offline_ns_parameters = np.zeros( 0 )
    M_ns_test = 0
    M_test_snapshots_matrix = np.zeros( 0 )
    M_test_parameters = np.zeros( 0 )

    M_N = 0
    M_basis = np.zeros( 0 )
    M_save_basis_functions = False
    M_save_file_basis_functions = "basis.txt"
    M_affineDecomposition = ad.AffineDecompositionHandler( )
    M_fem_problem = fm.fem_problem

    M_rb_handler = RbHandler( )

    def solve_reduced_problem( self, _param ):
        
        print( "Solving RB problem for parameter: " )
        print( _param )
        
        self.build_reduced_problem( _param )
        self.M_un = np.linalg.solve( self.M_An, self.M_fn )
        
        return
    
    def build_reduced_problem( self, _param ):
        
        N = self.M_N
        self.M_An = np.zeros( (N, N) )
        self.M_fn = np.zeros( N )
        self.M_un = np.zeros( N )
        
        for iQa in range( self.M_affineDecomposition.get_Qa( ) ):
            self.M_An = self.M_An + self.M_fem_problem.get_theta_a( _param, iQa ) * self.M_affineDecomposition.get_rb_affine_matrix( iQa )

        for iQf in range( self.M_affineDecomposition.get_Qf( ) ):
            self.M_fn = self.M_fn + self.M_fem_problem.get_theta_f( _param, iQf ) * self.M_affineDecomposition.get_rb_affine_vector( iQf )
        
        return
    
    def reconstruct_fem_solution( self, _un ):
    
        self.M_utildeh = np.zeros( (self.M_snapshots_matrix.shape[0], 1) )
        N = self.M_N
        
        assert _un.shape[0] == N
        
        self.M_utildeh = self.M_basis.dot( _un )
        
        return
        
    def print_rb_solution( self ):
        print( "\nThe RB solution is: " )
        print( self.M_un )
    
    def compute_rb_snapshots_error( self, _snapshot_number ):
        
        self.solve_reduced_problem( self.M_offline_ns_parameters[_snapshot_number, :] )
        self.reconstruct_fem_solution( self.M_un )
        
        error = self.M_utildeh
        error = error - self.M_snapshots_matrix[:, _snapshot_number]
        
        norm_of_error = np.linalg.norm( error ) / np.linalg.norm( self.M_snapshots_matrix[:, _snapshot_number] )
        
        print( "The norm of the error for snapshot %d is %g" % (_snapshot_number, norm_of_error) )
        
        return
    
    def compute_rb_test_snapshots_error( self, _snapshot_number ):
        
        self.solve_reduced_problem( self.M_test_parameters[_snapshot_number, :] )
        self.reconstruct_fem_solution( self.M_un )
        
        error = self.M_utildeh
        error = error - self.M_test_snapshots_matrix[:, _snapshot_number]
        
        norm_of_error = np.linalg.norm( error ) / np.linalg.norm( self.M_test_snapshots_matrix[:, _snapshot_number] )
        
        print( "The norm of the error for TEST snapshot %d is %g" % (_snapshot_number, norm_of_error) )
        
        return

    M_An = np.zeros( (0, 0) )
    M_fn = np.zeros( 0 )
    M_un = np.zeros( 0 )
    M_utildeh = np.zeros( 0 )






















