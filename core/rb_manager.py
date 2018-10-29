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
import fom_problem as fm
import parameter_handler as ph

class RbManager( ):
    
    def __init__( self, _affine_decomposition, _fom_problem, _parameter_handler ):
        
        self.set_affine_decomposition_handler( _affine_decomposition )
        self.set_fom_problem( _fom_problem )
        self.set_parameter_handler( _parameter_handler )
        
        return

    def import_snapshots_parameters( self, _input_file ):
        self.M_offline_ns_parameters = np.loadtxt( _input_file )
        return 

    def import_snapshots_matrix( self, _input_file ):
        self.M_snapshots_matrix = np.loadtxt( _input_file )
        self.M_ns = self.M_snapshots_matrix.shape[1]
        return 

    def import_snapshots( self, _input_file_snapshots, _input_file_snapshots_parameters ):

        self.import_snapshots_matrix( _input_file_snapshots )
        self.import_snapshots_parameters( _input_file_snapshots_parameters )
        
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
    
    def get_test_parameter( self, _iP ):
        return self.M_test_parameters[_iP, :]
    
    def get_snapshots_matrix( self, _fom_coordinates=np.array([]) ):
        if _fom_coordinates.shape[0]==0:
            return self.M_snapshots_matrix[:, :]
        else:
            return self.M_snapshots_matrix[_fom_coordinates.astype(int), :]
    
    def get_snapshot_function( self, _snapshot_number, _fom_coordinates=np.array([]) ):
        if self.M_get_test:
            return self.get_test_snapshot( _snapshot_number, _fom_coordinates )
        else:
            return self.get_snapshot( _snapshot_number, _fom_coordinates )   
        
    def get_snapshot( self, _snapshot_number, _fom_coordinates=np.array([]) ):
        
        if _fom_coordinates.shape[0]==0:
            return self.M_snapshots_matrix[:, _snapshot_number]
        else:
            return self.M_snapshots_matrix[_fom_coordinates.astype(int), _snapshot_number]
    
    def get_test_snapshot( self, _snapshot_number, _fom_coordinates=np.array([]) ):
    
        if _fom_coordinates.shape[0]==0:
            return self.M_test_snapshots_matrix[:, _snapshot_number]
        else:
            return self.M_test_snapshots_matrix[_fom_coordinates.astype(int), _snapshot_number]

    def get_number_of_snapshots( self ):
        return self.M_ns
    
    def get_basis( self, _fom_coordinates=np.array([]) ):
        if _fom_coordinates.shape[0]==0:
            return self.M_basis[:, :]
        else:
            return self.M_basis[_fom_coordinates.astype(int), :]
    
    def get_number_of_basis( self ):
        return self.M_N
    
    def get_fom_dimension( self ):
        return self.M_snapshots_matrix.shape[0]
    
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

    def set_fom_problem( self, _fom_problem ):
        self.M_fom_problem = _fom_problem
        return

    def set_parameter_handler( self, _parameter_handler ):
        self.M_parameter_handler = _parameter_handler
        return

    def reset_rb_approximation( self ):
        
        print( "Resetting RB approximation" )
        
        self.M_N = 0
        self.M_basis = np.zeros( 0 )
        self.M_affineDecomposition.reset_rb_approximation( )

    def build_rb_approximation( self, _ns, _tol = 10**(-5) ):
        
        self.reset_rb_approximation( )
        
        if self.M_ns < _ns :
            print( 'We miss some snalshots! I have only %d in memoory and I need to compute %d more.' % (self.M_ns, _ns-self.M_ns) )
            self.build_snapshots( _ns - self.M_ns )
        
        
        self.perform_pod( _tol )
        
        self.M_affineDecomposition.build_rb_affine_decompositions( self.M_basis )
        
        return

    M_verbose = False
    M_get_test = False 
    
    M_parameter_handler = ph.Parameter_handler( )
    
    M_ns = 0
    M_snapshots_matrix = np.zeros( ( 0, 0 ) )
    M_offline_ns_parameters = np.zeros( ( 0, 0 ) )
    M_ns_test = 0
    M_test_snapshots_matrix = np.zeros( 0 )
    M_test_parameters = np.zeros( 0 )

    M_N = 0
    M_basis = np.zeros( 0 )
    M_save_basis_functions = False
    M_save_file_basis_functions = "basis.txt"
    M_affineDecomposition = ad.AffineDecompositionHandler( )
    M_fom_problem = fm.fom_problem

    def get_rb_affine_matrix( self, _q ):
        return self.M_affineDecomposition.get_rb_affine_matrix( _q )

    def get_rb_affine_vector( self, _q ):
        return self.M_affineDecomposition.get_rb_affine_vector( _q )

    def solve_reduced_problem( self, _param ):
        
        if self.M_verbose == True:
            print( "Solving RB problem for parameter: " )
            print( _param )
        
        self.build_reduced_problem( _param )
        self.M_un = np.linalg.solve( self.M_An, self.M_fn )
        
        return self.M_un
    
    def build_reduced_problem( self, _param ):
        
        N = self.M_N
        self.M_An = np.zeros( (N, N) )
        self.M_fn = np.zeros( N )
        self.M_un = np.zeros( N )
        
        for iQa in range( self.M_affineDecomposition.get_Qa( ) ):
            self.M_An = self.M_An + self.M_fom_problem.get_theta_a( _param, iQa ) * self.get_rb_affine_matrix( iQa )

        for iQf in range( self.M_affineDecomposition.get_Qf( ) ):
            self.M_fn = self.M_fn + self.M_fom_problem.get_theta_f( _param, iQf ) * self.get_rb_affine_vector( iQf )
        
        return
    
    def reconstruct_fem_solution( self, _un ):
    
        self.M_utildeh = np.zeros( (self.M_snapshots_matrix.shape[0], 1) )
        
        assert _un.shape[0] == self.M_N
        
        self.M_utildeh = self.M_basis.dot( _un )
        
        return
        
    def get_utildeh( self ):
        return self.M_utildeh
    
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
        
#        norm_of_error = np.linalg.norm( error ) / np.linalg.norm( self.M_test_snapshots_matrix[:, _snapshot_number] )
        norm_of_error = np.sqrt( np.sum( error * error ) / \
            np.sum( self.M_test_snapshots_matrix[:, _snapshot_number] * self.M_test_snapshots_matrix[:, _snapshot_number] ) )
        
        print( "The norm of the error for TEST snapshot %d is %g" % (_snapshot_number, norm_of_error) )
        
        return norm_of_error

    M_An = np.zeros( (0, 0) )
    M_fn = np.zeros( 0 )
    M_un = np.zeros( 0 )
    M_utildeh = np.zeros( 0 )






















