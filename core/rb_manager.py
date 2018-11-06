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
        
        while cumulative_energy / total_energy < 1. - _tol**2:
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
    

class Mdeim( ):
    
    def __init__( self, _fom_problem ):
        self.M_fom_problem = _fom_problem
        return

    def reset_mdeim( self ):
        self.M_snapshots_matrix = np.zeros( ( 0, 0 ) )
        self.M_basis = np.zeros( ( 0, 0 ) )
        self.M_N = 0
        self.M_ns = 0
        
        return

    def build_mdeim_snapshots( self, _ns ):
        
        self.M_ns = _ns
        
        current_snapshots_number = self.M_snapshots_matrix.shape[1]
        
        for iNs in range( _ns ):
            
            self.M_fom_problem.generate_parameter( )
            param = self.M_fom_problem.get_parameter( )

            print('Choosing parameter %d for MDEIM' % iNs )

            AAA = self.M_fom_problem.assemble_fom_matrix( param )
            AA = np.array( AAA['A'] )
 
            if current_snapshots_number == 0:
                self.M_snapshots_matrix = np.zeros( ( AA.shape[0], _ns ) )
                self.M_row_map = AA[:, 0]
                self.M_col_map = AA[:, 1]

            self.M_snapshots_matrix[:, iNs] = AA[:, 2]
            current_snapshots_number = current_snapshots_number + 1
       
        return

    def build_mdeim_basis( self, _ns, _tol ):
        
        self.reset_mdeim( )
        
        self.build_mdeim_snapshots( _ns )
        
        pod = ProperOrthogonalDecompostion( )
        
        pod.perform_pod( self.M_snapshots_matrix, _tol )
        
        self.M_basis = pod.get_basis( )
    
        self.M_N = self.M_basis.shape[1]

        print( 'MDEIM BASIS' )
        print( self.M_basis )
 
        if self.M_save_mdeim_basis == True:
            output_file = open( self.M_save_file_basis_functions + '_' + '20' + '.txt', 'w' )

            for iV in range( self.M_basis.shape[0] ):
                for iB in range( self.M_basis.shape[1] ):
                    output_file.write( "%.10g" % self.M_basis[iV, iB] )

                    if iB < self.M_basis.shape[1] - 1:
                        output_file.write( " " )
                    else:
                        output_file.write( "\n" )

            output_file.close( )


            for iB in range( self.M_basis.shape[1] ):
                output_file = open( self.M_save_file_basis_functions + '_A' + '20' + '_' + str(iB) + '.txt', 'w' )
                
                for iV in range( self.M_basis.shape[0] ):
                    output_file.write( "%d %d %.10g\n" % (self.M_row_map[iV], self.M_col_map[iV], self.M_basis[iV, iB]) )

                output_file.close( )

    def get_num_mdeim_basis( self ):
        return self.M_N

    M_fom_problem = 0
    M_snapshots_matrix = np.zeros( ( 0, 0 ) )
    M_basis = np.zeros( ( 0, 0 ) )
    M_N = 0
    M_ns = 0

    M_row_map = np.zeros( ( 0, 0 ) )
    M_col_map = np.zeros( ( 0, 0 ) )

    M_save_file_basis_functions = "mdeim_basis"
    M_save_mdeim_basis = True


class RbManager( ):
    
    def __init__( self, _affine_decomposition, _fom_problem ):
        
        self.set_affine_decomposition_handler( _affine_decomposition )
        self.set_fom_problem( _fom_problem )
        
        return

    def import_snapshots_parameters( self, _input_file ):
        self.M_offline_ns_parameters = np.loadtxt( _input_file )
        return 

    def import_snapshots_matrix( self, _input_file ):
        self.M_snapshots_matrix = np.loadtxt( _input_file )
        self.M_ns = self.M_snapshots_matrix.shape[1]
        return 

    def import_basis_matrix( self, _input_file ):
        self.M_basis = np.loadtxt( _input_file )
        self.M_N = self.M_basis.shape[1]
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
   
    # _ns is the number of snapshots to be added to the snapshots matrix
    def build_snapshots( self, _new_snapshots ):
        
        current_snapshots_number = self.M_snapshots_matrix.shape[1]

        if current_snapshots_number == 0:
            print( "There are no snapshots at the moment, therefore we need to create everything from scratch" )

        num_parameters = self.M_fom_problem.get_num_parameters( )
        new_parameters = np.zeros( ( _new_snapshots, num_parameters ) )
        
        for iS in range( _new_snapshots ):
            
            self.M_fom_problem.generate_parameter( )
            new_parameters[iS, :] = self.M_fom_problem.get_parameter( )
            
            print( "Considering the parameter %d " % iS )
            print( new_parameters[iS, :] )

            sol = self.M_fom_problem.solve_fom_problem( new_parameters[iS, :] )
            
            if current_snapshots_number == 0:
                self.M_snapshots_matrix = np.zeros( ( len(sol['u']), _new_snapshots ) )

            u = np.array( sol['u'] )

            self.M_snapshots_matrix[:, iS] = u[:, 0]
            current_snapshots_number = current_snapshots_number + 1
        
        self.M_ns = self.M_ns + _new_snapshots
        
        
        if self.M_save_offline_structures == True:
            output_file = open( self.M_save_file_snapshots_functions, 'w+' )
                
            for iNs in range( self.M_snapshots_matrix.shape[0] ):
                for iNh in range( self.M_snapshots_matrix.shape[1] ):
                    output_file.write( "%.10g" % self.M_snapshots_matrix[iNs, iNh] )
    
                    if iNh < self.M_snapshots_matrix.shape[1] - 1:
                        output_file.write( " " )
                    else:
                        output_file.write( "\n" )
            
            output_file.close( )
            
            output_file = open( self.M_save_file_offline_parameters, 'w+' )

            for iNs in range( new_parameters.shape[0] ):
                for iNh in range( new_parameters.shape[1] ):
                    output_file.write( "%.10g" % new_parameters[iNs, iNh] )
    
                    if iNh < new_parameters.shape[1] - 1:
                        output_file.write( " " )
                    else:
                        output_file.write( "\n" )
            
            output_file.close( )

        return

    def perform_pod( self, _tol = 10**(-5) ):

        pod = ProperOrthogonalDecompostion( )
        
        pod.perform_pod( self.M_snapshots_matrix, _tol )
        
        self.M_basis = pod.get_basis( )
    
        self.M_N = self.M_basis.shape[1]
    
        if self.M_save_offline_structures == True:
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

    def reset_rb_approximation( self ):
        
        print( "Resetting RB approximation" )
        
        self.M_N = 0
        self.M_basis = np.zeros( 0 )
        self.M_affineDecomposition.reset_rb_approximation( )

    def build_rb_approximation( self, _ns, _tol = 10**(-5) ):
        
        self.reset_rb_approximation( )
        
        if self.M_ns < _ns :
            print( 'We miss some snalshots! I have only %d in memory and I need to compute %d more.' % (self.M_ns, _ns-self.M_ns) )
            self.build_snapshots( _ns - self.M_ns )
        
        self.perform_pod( _tol )
        
        self.M_affineDecomposition.build_rb_affine_decompositions( self.M_basis, self.M_fom_problem )
        
        if self.M_save_offline_structures == True:
            self.M_affineDecomposition.save_rb_affine_decomposition( self.M_save_file_affine_components )

        return

    # to save the offline structures in the given locations
    def save_offline_structures( self, _snapshot_matrix, _basis_matrix, _affine_components, _offline_parameters ):
        
        self.M_save_offline_structures       = True
        self.M_save_file_basis_functions     = _basis_matrix
        self.M_save_file_snapshots_functions = _snapshot_matrix
        self.M_save_file_affine_components   = _affine_components
        self.M_save_file_offline_parameters  = _offline_parameters

        return
    
    def import_offline_structures( self, _snapshot_matrix, _basis_matrix, _affine_components ):
        
        self.import_snapshots_matrix( _snapshot_matrix )
        self.import_basis_matrix( _basis_matrix )
        self.M_affineDecomposition.import_affine_components( _affine_components )
        
        return
    
    M_verbose = False
    M_get_test = False 
    
    M_ns = 0
    M_snapshots_matrix = np.zeros( ( 0, 0 ) )
    M_offline_ns_parameters = np.zeros( ( 0, 0 ) )
    M_ns_test = 0
    M_test_snapshots_matrix = np.zeros( 0 )
    M_test_parameters = np.zeros( 0 )

    M_N = 0
    M_basis = np.zeros( 0 )
    M_affineDecomposition = ad.AffineDecompositionHandler( )
    M_fom_problem = fm.fom_problem

    M_save_offline_structures = False
    M_save_file_basis_functions = "basis.txt"
    M_save_file_offline_parameters = "offline_parameters.data"

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

    def test_rb_solver( self, _n_test ):
        
        all_errors = 0.;
        
        for iP in range( _n_test ):
            
            new_param = self.M_fom_problem.generate_parameter( )
            new_param = self.M_fom_problem.get_parameter( )

            print( "New parameter %d " %( iP ) )
            print( new_param )

            self.solve_reduced_problem( new_param )
            self.reconstruct_fem_solution( self.M_un )
            
            uhh = self.M_fom_problem.solve_fom_problem( new_param )
            uh = np.array( uhh['u'] )
            uh = uh[:, 0]
            error = self.M_utildeh
            error = error - uh
            
            norm_of_error = np.sqrt( np.sum( error * error ) / np.sum( uh * uh )  )
            all_errors = all_errors + norm_of_error
            
            print( "The error is %e \n\n" % norm_of_error )

        avg_error = all_errors / _n_test

        print( "The average error is %E" % avg_error )
        
        return avg_error


    M_An = np.zeros( (0, 0) )
    M_fn = np.zeros( 0 )
    M_un = np.zeros( 0 )
    M_utildeh = np.zeros( 0 )

















