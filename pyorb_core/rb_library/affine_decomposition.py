#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 14:41:45 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np
import pyorb_core.algebraic_utils as alg_ut

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

            #  if the matrix indices start from 1, we rescale them
            if np.min( self.M_feAffineAq[iQa][:, 0:2] ) > 0 :
                self.M_feAffineAq[iQa][:, 0:2] = self.M_feAffineAq[iQa][:, 0:2] - 1

        return

    # _input_file should be the string that have in common the affine matrices
    def import_affine_vectors( self, _input_file ):

        Qf = self.M_affineDecomposition.get_Qf( )

        assert Qf > 0

        self.M_feAffineFq = []

        for iQf in range( Qf ):
            self.M_feAffineFq.append( np.loadtxt( _input_file + str(iQf) + '.txt' ) )   # importing vectors

        return

    def import_rb_affine_matrices( self, _input_file ):

        Qa = self.M_affineDecomposition.get_Qa( )

        assert Qa > 0

        self.M_rbAffineAq = []

        for iQa in range( Qa ):
            self.M_rbAffineAq.append( np.loadtxt( _input_file + str(iQa) + '.txt' ) )   # importing rb matrix

        return

    # _input_file should be the string that have in common the affine matrices
    def import_rb_affine_vectors( self, _input_file ):

        Qf = self.M_affineDecomposition.get_Qf( )

        assert Qf > 0

        self.M_rbAffineFq = []

        for iQf in range( Qf ):
            self.M_rbAffineFq.append( np.loadtxt( _input_file + str(iQf) + '.txt' ) )   # importing rb vectors

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


    def set_affine_a( self, _feAffineAq ):
        self.M_feAffineAq = _feAffineAq


    def set_affine_f( self, _feAffineFq ):
        self.M_feAffineFq = _feAffineFq


    def build_rb_affine_decompositions( self, _basis, _fom_problem, _build_rb_tpl=False ):

        N = _basis.shape[1]

        Qf = self.get_Qf( )
        Qa = self.get_Qa( )

        if _build_rb_tpl is not True:
          if self.check_set_fom_arrays( ) == False:

              print( "Importing FOM affine arrays " )

              if len( self.M_feAffineFq ) < Qf:
                  print( "I am importing f affine arrays " )

                  fff = _fom_problem.retrieve_fom_affine_components( 'f', Qf - len( self.M_feAffineFq ) )
                  starting_Qf = len( self.M_feAffineFq )

                  for iQf in range( Qf - starting_Qf ):
                    self.M_feAffineFq.append( np.array( fff['f' + str(iQf)] ) )

              if len( self.M_feAffineAq ) < Qa:
                  print( "I am importing A affine arrays starting from %d and to %d " % (len( self.M_feAffineAq ), Qa) )

                  AAA = _fom_problem.retrieve_fom_affine_components( 'A', Qa - len( self.M_feAffineAq ) )
                  starting_Qa = len( self.M_feAffineAq )

                  for iQa in range( Qa - starting_Qa ):
                      print( "I am importing A affine array %d " % (iQa + starting_Qa) )

                      self.M_feAffineAq.append( AAA['A' + str(iQa)] )

              if len( self.M_feAffineAq ) < Qa:
                  print( "I am importing A - Jacobian affine arrays starting from %d and to %d " % (len( self.M_feAffineAq ), Qa) )

                  AAA = _fom_problem.retrieve_fom_affine_components( 'Aj', Qa - len( self.M_feAffineAq ) )
                  starting_Qa = len( self.M_feAffineAq )

                  for iQa in range( Qa - starting_Qa ):
                      print( "I am importing A affine array %d " % (iQa + starting_Qa) )

                      self.M_feAffineAq.append( AAA['A' + str(iQa)] )

          else:
              print( "Already set the FOM affine arrays " )

          for iQf in range( Qf ):
            self.M_rbAffineFq.append( np.zeros( N ) )
            self.M_rbAffineFq[iQf] = _basis.T.dot( self.M_feAffineFq[iQf] )

          for iQa in range( Qa ):
            Av = alg_ut.sparse_matrix_vector_mul( self.M_feAffineAq[iQa], _basis )
            self.M_rbAffineAq.append( np.zeros( (N, N) ) )
            self.M_rbAffineAq[iQa] = _basis.T.dot( Av )

        elif _build_rb_tpl is True:
 
          print( "Importing directly the RB arrays from TPL " )        
  
          rbAffineFq_components = _fom_problem.retrieve_rb_affine_components( 'f' )
          
          for iQf in range( len( rbAffineFq_components ) ):
            self.M_rbAffineFq.append( rbAffineFq_components[iQf] ) 

          rbAffineAq_components = _fom_problem.retrieve_rb_affine_components( 'A' )

          for iQa in range( len( rbAffineAq_components ) ):
            self.M_rbAffineAq.append( rbAffineAq_components[iQa] )

          if len( self.M_feAffineAq ) < Qa:
            rbAffineAjq_components = _fom_problem.retrieve_rb_affine_components( 'Aj' )

            for iQaj in range( len( rbAffineAjq_components ) ):
              self.M_rbAffineAq.append( rbAffineAjq_components[iQaj] )

        print( 'Finished to build the RB affine arrays' )

        return

    def check_set_fom_arrays( self ):
        return len( self.M_feAffineAq ) >= self.get_Qa( ) and len( self.M_feAffineFq ) >= self.get_Qf( )

    def save_rb_affine_decomposition( self, _file_name ):

        Qf = self.get_Qf( )
        Qa = self.get_Qa( )

        for iQa in range( Qa ):
            output_file = open( _file_name + '_A' + str( iQa ) + '.txt', 'w+' )

            for iN in range( self.M_rbAffineAq[iQa].shape[0] ):
                for jN in range( self.M_rbAffineAq[iQa].shape[1] ):
                    output_file.write( "%.10g" % self.M_rbAffineAq[iQa][iN, jN] )

                    if jN < self.M_rbAffineAq[iQa].shape[1] - 1:
                        output_file.write( " " % self.M_rbAffineAq[iQa][iN, jN] )
                    else:
                        output_file.write( "\n" % self.M_rbAffineAq[iQa][iN, jN] )

            output_file.close( )

        for iQf in range( Qf ):
            output_file = open( _file_name + '_f' + str( iQf ) + '.txt', 'w+' )

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
            self.M_rbAffineAq.append( np.loadtxt( _affine_components + '_A' + str( iQa ) + '.txt' ) )

        for iQf in range( Qf ):
            self.M_rbAffineFq.append( np.loadtxt( _affine_components + '_f' + str( iQf ) + '.txt' ) )

        return

    M_feAffineAq = []
    M_feAffineFq = []

    M_rbAffineAq = []
    M_rbAffineFq = []

    M_affineDecomposition = AffineDecomposition( )
