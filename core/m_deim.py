#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  8 18:21:34 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np
import proper_orthogonal_decomposition as podec

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
        
        pod = podec.ProperOrthogonalDecompostion( )
        
        pod.perform_pod( self.M_snapshots_matrix, _tol )
        
        self.M_basis = pod.get_basis( )
    
        self.M_N = self.M_basis.shape[1]

        print( 'MDEIM BASIS' )
        print( self.M_basis )
 
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
