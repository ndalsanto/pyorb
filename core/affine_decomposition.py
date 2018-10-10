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

    # _input_file should be the string that have in common the affine matrices
    def import_affine_matrices( self, _input_file ):
        
        Qa = self.M_affineDecomposition.get_Qa( )

        assert Qa > 0
        
        for iQa in range( Qa ):
            self.M_feAffineAq.append( np.loadtxt( _input_file + str(iQa) + '.txt' ) )   # importing matrix in sparse format
        
        return
        
    # _input_file should be the string that have in common the affine matrices
    def import_affine_vectors( self, _input_file ):
        
        Qf = self.M_affineDecomposition.get_Qf( )

        assert Qf > 0
        
        for iQf in range( Qf ):
            self.M_feAffineFq.append( np.loadtxt( _input_file + str(iQf) + '.txt' ) )   # importing vectors
        
        return
        
#    def resize_rb_arrays( self, _N ):
#        
#        for iQf in range( self.M_qf ):
#            self.M_rbAffineFq.append(  )
#        
#        for iQa in range( self.M_qa ):
#            self.M_rbAffineAq.append( np.zeros( (_N, _N) ) )   # importing matrix in sparse format
#        
#        return
    
    def print_ad_summary( self ):
        
        self.M_affineDecomposition.print_ad_summary( )
        
        return 
    
    def build_rb_affine_decompositions( self, _basis ):
        
        N = _basis.shape[1]
        
        Qf = self.get_Qf( )
        
        for iQf in range( Qf ):
            
            self.M_rbAffineFq.append( np.zeros( N ) )
            self.M_rbAffineFq[iQf] = _basis.T.dot( self.M_feAffineFq[iQf] )
        
        Qa = self.get_Qa( )

        for iQa in range( Qa ):
            Av = sparse_matrix_vector_mul( self.M_feAffineAq[iQa], _basis )
            self.M_rbAffineAq.append( np.zeros( (N, N) ) )
            self.M_rbAffineAq[iQa] = _basis.T.dot( Av )
            
            print( "\n\n Affine component %d " % iQa )
            print( self.M_rbAffineAq[iQa] )
        
        return

    M_feAffineAq = []
    M_feAffineFq = []

    M_rbAffineAq = []
    M_rbAffineFq = []
    
    M_affineDecomposition = AffineDecomposition( )



















