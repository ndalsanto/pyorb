#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 14:41:45 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np

# mat is a matrix in COO format, vec is a vector, such that the results is Av = mat * vec

def sparse_matrix_vector_mul( mat, vec ):

    Av = np.zeros( vec.shape  )

    nnz = mat.shape[0]
    
    for i in range( nnz ):
        Av[ int(mat[i, 0])-1, :] = Av[ int(mat[i, 0])-1, :] + mat[i, 2] * vec[ int(mat[i, 1])-1, : ]
    
    return Av


class AffineDecomposition( ):

    def __init__( self ):
        return
    
    def get_Qa( self ):
        return self.M_qa

    def get_Qf( self ):
        return self.M_qf

    def set_Q( self, _qa, _qf ):
        self.M_qa = _qa
        self.M_qf = _qf
        return

    def print_ad_summary( self ):
        
        print( "Number of affine decomposition matrices %d" % self.M_qa )
        print( "Number of affine decomposition vectors  %d"  % self.M_qf )
        
        return 
    
    M_qa = 0
    M_qf = 0


class AffineDecompositionHandler( ):
    
    def __init__( self ):
        return
    
    def get_Qa( self ):
        return self.M_affineDecomposition.get_Qa( )

    def get_Qf( self ):
        return self.M_affineDecomposition.get_Qf( )

    def set_Q( self, _qa, _qf ):
        self.M_affineDecomposition.set_Q( _qa, _qf )
        return
 
    def get_affine_matrix( self, _q ):
        return self.M_feAffineAq[_q]

    def get_affine_vector( self, _q ):
        return self.M_feAffineFq[_q]
    
    def get_rb_affine_matrix( self, _q ):
        return self.M_rbAffineAq[_q]

    def get_rb_affine_vector( self, _q ):
        return self.M_rbAffineFq[_q]

    # _input_file should be the string that have in common the affine matrices
    def import_affine_matrices( self, _input_file ):
        
        Qa = self.M_affineDecomposition.get_Qa( )

        assert Qa > 0
        
        self.M_feAffineAq = []
        
        for iQa in range( Qa ):
            self.M_feAffineAq.append( np.loadtxt( _input_file + str(iQa) + '.txt' ) )   # importing matrix in sparse format
        
        return
        
    # _input_file should be the string that have in common the affine matrices
    def import_affine_vectors( self, _input_file ):
        
        Qf = self.M_affineDecomposition.get_Qf( )

        assert Qf > 0

        self.M_feAffineFq = []
        
        for iQf in range( Qf ):
            self.M_feAffineFq.append( np.loadtxt( _input_file + str(iQf) + '.txt' ) )   # importing vectors
        
        return
        
    def print_ad_summary( self ):
        
        self.M_affineDecomposition.print_ad_summary( )
        
        return 
    
    def print_affine_components( self ):
    
        Qf = self.get_Qf( )
        Qa = self.get_Qa( )

        for iQf in range( Qf ):
            print( '\nRB rhs affine components %d \n' % iQf )
            print( self.M_rbAffineFq[iQf] ) 
            
        for iQa in range( Qa ):
            print( '\nRB mat affine components %d \n' % iQa )
            print( self.M_rbAffineAq[iQa] ) 
                    
        return 

    def reset_rb_approximation( self ):
        self.M_rbAffineFq = []
        self.M_rbAffineAq = []

    
    def build_rb_affine_decompositions( self, _basis, _fom_problem ):
        
        N = _basis.shape[1]
        
        Qf = self.get_Qf( )
        Qa = self.get_Qa( )

        if self.check_set_fom_arrays( ) == False:

            print( "Importing FOM affine arrays " )

            fff = _fom_problem.retrieve_fe_affine_components( 'f' )

            for iQf in range( Qf ):
                self.M_feAffineFq.append( np.array( fff['f' + str(iQf)] ) )
                self.M_feAffineFq[iQf] = self.M_feAffineFq[iQf][:, 0]
                
            AAA = _fom_problem.retrieve_fe_affine_components( 'A' )

            for iQa in range( Qa ):
                self.M_feAffineAq.append( np.array( AAA['A' + str(iQa)] ) )

        for iQf in range( Qf ):
            self.M_rbAffineFq.append( np.zeros( N ) )
            self.M_rbAffineFq[iQf] = _basis.T.dot( self.M_feAffineFq[iQf] )
            
        for iQa in range( Qa ):
            Av = sparse_matrix_vector_mul( self.M_feAffineAq[iQa], _basis )
            self.M_rbAffineAq.append( np.zeros( (N, N) ) )
            self.M_rbAffineAq[iQa] = _basis.T.dot( Av )
        
        return


    def check_set_fom_arrays( self ):
        return len( self.M_feAffineAq ) > 0 and len( self.M_feAffineFq ) > 0

    def save_rb_affine_decomposition( self, _file_name ):
        
        Qf = self.get_Qf( )
        Qa = self.get_Qa( )

        for iQa in range( Qa ):
            output_file = open( _file_name + 'A' + str( iQa ), 'w+' )
                    
            for iN in range( self.M_rbAffineAq[iQa].shape[0] ):
                for jN in range( self.M_rbAffineAq[iQa].shape[1] ):
                    output_file.write( "%.10g" % self.M_rbAffineAq[iQa][iN, jN] )
    
                    if jN < self.M_rbAffineAq[iQa].shape[1] - 1:
                        output_file.write( " " % self.M_rbAffineAq[iQa][iN, jN] )
                    else:
                        output_file.write( "\n" % self.M_rbAffineAq[iQa][iN, jN] )
            
            output_file.close( )
    
        for iQf in range( Qf ):
            output_file = open( _file_name + '_f' + str( iQf ), 'w+' )
                    
            for iN in range( self.M_rbAffineFq[iQf].shape[0] ):
                output_file.write( "%.10g" % self.M_rbAffineFq[iQf][iN] )
                output_file.write( " " % self.M_rbAffineFq[iQf][iN] )
                output_file.write( "\n" % self.M_rbAffineFq[iQf][iN] )
            
            output_file.close( )
        
        return

    def import_affine_components( self, _affine_components ):
        
        self.M_rbAffineAq = []
        self.M_rbAffineFq = []

        Qf = self.get_Qf( )
        Qa = self.get_Qa( )

        for iQa in range( Qa ):
            self.M_rbAffineAq.append( np.loadtxt( _affine_components + 'A' + str( iQa ) ) )
            
        for iQf in range( Qf ):
            self.M_rbAffineFq.append( np.loadtxt( _affine_components + '_f' + str( iQf ) ) )
            
        return

    M_feAffineAq = []
    M_feAffineFq = []

    M_rbAffineAq = []
    M_rbAffineFq = []
    
    M_affineDecomposition = AffineDecomposition( )















