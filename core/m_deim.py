#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 18:21:34 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np
import proper_orthogonal_decomposition as podec
import random

class Deim( ):
    
    def __init__( self, _fom_problem ):
        self.M_fom_problem = _fom_problem
        self.M_current_param = np.zeros( _fom_problem.get_num_parameters( ) )
        return

    def reset_deim( self ):
        self.M_snapshots_matrix = np.zeros( ( 0, 0 ) )
        self.M_basis = np.zeros( ( 0, 0 ) )
        self.M_N = 0
        self.M_ns = 0
        
        return

    def build_deim_snapshots( self, _ns ):
        
        print('You need to code this !! ' )

#        self.M_ns = _ns
        
#        current_snapshots_number = self.M_snapshots_matrix.shape[1]
        
#        for iNs in range( _ns ):
#            
#            self.M_fom_problem.generate_parameter( )
#            param = self.M_fom_problem.get_parameter( )
#
#            print('Choosing parameter %d for MDEIM' % iNs )
#
#            AAA = self.M_fom_problem.assemble_fom_matrix( param )
#            AA = np.array( AAA['A'] )
# 
#            if current_snapshots_number == 0:
#                self.M_snapshots_matrix = np.zeros( ( AA.shape[0], _ns ) )
#                self.M_row_map = AA[:, 0]
#                self.M_col_map = AA[:, 1]
#
#            self.M_snapshots_matrix[:, iNs] = AA[:, 2]
#            current_snapshots_number = current_snapshots_number + 1
#       
        return

    def build_deim_basis( self, _ns, _tol ):
        
        pod = podec.ProperOrthogonalDecompostion( )
        
        pod.perform_pod( self.M_snapshots_matrix, _tol )
        
        self.M_basis = pod.get_basis( )
    
        self.M_N = self.M_basis.shape[1]

#        print( 'MDEIM BASIS' )
#        print( self.M_basis )
 
        if self.M_save_mdeim_basis == True:
            output_file = open( self.M_save_file_basis_functions + '_' + str(self.M_fom_problem.M_fom_specifics['number_of_elements']) + '.txt', 'w' )

            for iV in range( self.M_basis.shape[0] ):
                for iB in range( self.M_basis.shape[1] ):
                    output_file.write( "%.10g" % self.M_basis[iV, iB] )

                    if iB < self.M_basis.shape[1] - 1:
                        output_file.write( " " )
                    else:
                        output_file.write( "\n" )

            output_file.close( )

            for iB in range( self.M_basis.shape[1] ):
                output_file = open( self.M_save_file_basis_functions + '_A' + str(self.M_fom_problem.M_fom_specifics['number_of_elements']) + '_' + str(iB) + '.txt', 'w' )
                
                for iV in range( self.M_basis.shape[0] ):
                    output_file.write( "%d %d %.10g\n" % (self.M_row_map[iV], self.M_col_map[iV], self.M_basis[iV, iB]) )

                output_file.close( )

    def perform_deim( self, _ns, _tol ):
        
        self.reset_deim( )
        
        self.build_deim_snapshots( _ns )

        self.build_deim_basis( self, _ns, _tol )
        
        return

    def get_num_basis( self ):
        return self.M_N

    def print_reduced_indices( self ):
        print( self.M_reduced_indices )
        return

    def identify_reduced_indeces( self ):
        
        self.M_reduced_indices = np.zeros( self.M_N )
        
        self.M_reduced_indices[0] = np.argmax( np.abs( self.M_basis[:, 0] ) )
        
        res = self.M_basis[ self.M_reduced_indices[0].astype(int), 1] \
            / self.M_basis[ self.M_reduced_indices[0].astype(int), 0 ]
        
        r = self.M_basis[:, 1] - self.M_basis[:, 0].dot( res )

        self.M_reduced_indices[1] = np.argmax( abs( r ) ).astype( int )

        for iB in range(2, self.M_N):
            
            res = np.linalg.solve( self.M_basis[ self.M_reduced_indices[0:iB].astype(int), 0:iB ], 
                                   self.M_basis[ self.M_reduced_indices[0:iB].astype(int), iB])
            
            r = self.M_basis[:, iB] - self.M_basis[:, 0:iB].dot( res.T )
            
            self.M_reduced_indices[iB] = np.argmax( np.abs( r ) )
    
        self.M_reduced_indices = self.M_reduced_indices.astype( int )
    
        self.M_interpolation_matrix = self.M_basis[ self.M_reduced_indices[0:self.M_N].astype(int), :]
    
        return
    
    def get_basis( self ):
        
        return self.M_basis

    M_fom_problem = 0
    M_snapshots_matrix = np.zeros( ( 0, 0 ) )
    M_basis = np.zeros( ( 0, 0 ) )
    M_N = 0
    M_ns = 0
    M_offline_parameters = np.zeros( ( 0, 0 ) )

    M_interpolation_matrix = np.zeros( ( 0, 0 ) )
    M_reduced_indices = np.zeros( (0, 0) )

    M_save_file_basis_functions = "deim_basis"
    M_save_deim_basis = True
    
    M_current_theta = np.zeros( 0 )
    M_current_param = np.zeros( 0 )


class Mdeim( Deim ):
    
    def __init__( self, _fom_problem ):
        Deim.__init__( self, _fom_problem )
        self.M_save_file_basis_functions = "mdeim_basis"
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

        self.M_fom_problem.generate_parameter( )

        self.M_offline_parameters = np.zeros( (_ns, self.M_fom_problem.get_num_parameters( ) ) )
        
        for iNs in range( _ns ):
            
            random.seed( (iNs + 13) * 1013 )
            
            self.M_fom_problem.generate_parameter( )
            param = self.M_fom_problem.get_parameter( )

            self.M_offline_parameters[iNs, :] = param

            np.set_printoptions(precision=12)
            print( param )

            AAA = self.M_fom_problem.assemble_fom_matrix( param )
            AA = np.array( AAA['A'] )
 
            if current_snapshots_number == 0:
                self.M_snapshots_matrix = np.zeros( ( AA.shape[0], _ns ) )
                self.M_row_map = AA[:, 0].astype( int )
                self.M_col_map = AA[:, 1].astype( int )

            self.M_snapshots_matrix[:, iNs] = AA[:, 2]
            current_snapshots_number = current_snapshots_number + 1
       
        return

    def identify_reduced_mat_indeces( self ):
        
        self.M_reduced_indices_mat = np.zeros( (self.M_N, 2) )
        
        self.M_reduced_indices_mat[:, 0] = self.M_row_map[ self.M_reduced_indices ]
        self.M_reduced_indices_mat[:, 1] = self.M_col_map[ self.M_reduced_indices ]
        
        self.M_reduced_indices_mat = self.M_reduced_indices_mat.astype( int )
        
        return

    def identify_reduced_elements( self ):
        
        self.M_reduced_elements = self.M_fom_problem.find_mdeim_elements_fom_specifics( self.M_reduced_indices_mat )

        return

    def print_reduced_indices_mat( self ):
        
        print( self.M_reduced_indices_mat )
        
        return

    def print_reduced_elements( self ):
        
        print( self.M_reduced_elements )
        
        return

    def build_mdeim_basis( self, _ns, _tol ):
        return

    def perform_mdeim( self, _ns, _tol ):
        
        self.reset_mdeim( )
        
        self.build_mdeim_snapshots( _ns )
        
        self.compute_offline_mdeim( _ns, _tol )
        
        return

    def compute_offline_mdeim( self, _ns, _tol ):
        
        self.build_deim_basis( _ns, _tol )

        self.identify_reduced_indeces( )

        self.identify_reduced_mat_indeces( )

        self.identify_reduced_elements( )
        
        return

    def get_num_mdeim_basis( self ):
        return self.M_N

    def compute_theta_coefficients( self, _param ):
        
        rhs = self.M_fom_problem.assemble_fom_matrix( _param, _elements=self.M_reduced_elements, \
                                                      _indices=self.M_reduced_indices_mat )
        self.M_current_theta = np.linalg.solve( self.M_interpolation_matrix, rhs[:, 2] )
        
        return self.M_current_theta
    
    def compute_theta_coefficients_q( self, _param, _q ):
        
        if (self.M_current_param != _param).all():
            print( 'Recomputing for new paramter!' )
            self.M_current_param = _param + np.zeros( _param.shape )
            self.compute_theta_coefficients( _param )
        
        return self.M_current_theta[_q]
        
    
    def compute_theta_bounds( self, _n_tests ):
        
        self.M_fom_problem.generate_parameter( )
        param = self.M_fom_problem.get_parameter( )

        min_thetas = self.compute_theta_coefficients( param )
        max_thetas = min_thetas
        
        for iT in range( 1, _n_tests ):
            self.M_fom_problem.generate_parameter( )
            param = self.M_fom_problem.get_parameter( )
            this_theta = self.compute_theta_coefficients( param )

            max_thetas = np.maximum( this_theta, max_thetas )
            min_thetas = np.minimum( this_theta, min_thetas )

        return min_thetas, max_thetas
    
    def get_basis_list( self ):
        l = []
        
        for iB in range( self.M_N ):
            
            basis = np.zeros( (self.M_basis.shape[0], 3) )
            basis[:, 0] = self.M_row_map
            basis[:, 1] = self.M_col_map
            basis[:, 2] = self.M_basis[:, iB]

            l.append( basis )
            
        return l

    M_row_map = np.zeros( ( 0, 0 ) )
    M_col_map = np.zeros( ( 0, 0 ) )

    M_reduced_indices_mat = np.zeros( ( 0, 0 ) )
    M_reduced_elements = np.zeros( ( 0, 0 ) )

    M_save_mdeim_basis = False






