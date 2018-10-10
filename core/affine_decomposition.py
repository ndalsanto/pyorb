#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 14:41:45 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np

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
        return self.M_qa

    def get_Qf( self ):
        return self.M_qf

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
        
    def resize_rb_arrays( self, _N ):
        
        for iQf in range( self.M_qf ):
            self.M_rbAffineFq.append( np.zeros( _N ) )
        
        for iQa in range( self.M_qa ):
            self.M_rbAffineAq.append( np.zeros( (_N, _N) ) )   # importing matrix in sparse format
        
        return
    
    def print_ad_summary( self ):
        
        self.M_affineDecomposition.print_ad_summary( )
        
        return 
    
    def build_rb_affine_decompositions( self ):
        
        
        
        return

    M_feAffineAq = []
    M_feAffineFq = []

    M_rbAffineAq = []
    M_rbAffineFq = []
    
    M_affineDecomposition = AffineDecomposition( )



















