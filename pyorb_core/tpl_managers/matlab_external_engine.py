#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 13:51:49 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import matlab.engine
import numpy as np
from pyorb_core.tpl_managers import external_engine as ee

class matlab_external_engine( ee.external_engine ):
    
    def __init__( self, _engine_type, _library_path ):
        
        ee.external_engine.__init__( self, _engine_type, _library_path )
        
        return

    def start_specific_engine( self ):
        
        self.M_engine = matlab.engine.start_matlab( )
        self.M_engine.addpath( self.M_engine.genpath( self.M_library_path ) )

        print( 'Successfully started matlab engine and corresponding FOM library %s ' % self.M_library_path )

        return

    def quit_specific_engine( self ):

        self.M_engine.quit( )

        print( 'Successfully quitted matlab engine' )

        return

    def convert_parameter( self, _param ):

        return self.convert_double( _param )

    def convert_double( self, _np_array ):

        return matlab.double( _np_array.tolist() )

    def convert_indices( self, _indices ):

        return matlab.int64(_indices.tolist())

    def solve_parameter( self, _param, _fom_specifics ):

        u = self.M_engine.solve_parameter( self.convert_parameter( _param ), _fom_specifics )

        sol = np.array( u['u'] )

        return sol[:, 0]

    # normally a MATLAB application can directly provide a dictionary with all the affine components 
    def build_fom_affine_components( self, _operator, _num_affine_components, _fom_specifics ):

        print( 'Building affine components for operator %c' % _operator )
        
        affine_components = self.M_engine.build_fom_affine_components( _operator, _fom_specifics )

        # rescale the matrix indices so that the counting starts from 0 (and not from 1 as in MATLAB)
        if _operator == 'A':
            matrix_affine = { }
            for iQa in range( _num_affine_components ):
                matrix_affine['A' + str(iQa)] = np.array( affine_components['A' + str(iQa)] )
                matrix_affine['A' + str(iQa)][:, 0:2] = matrix_affine['A' + str(iQa)][:, 0:2] - 1
                
            affine_components = matrix_affine

        # resetting to just one-size array
        if _operator == 'f':
            rhs_affine = {}
            for iQf in range( _num_affine_components ):
                rhs_affine['f' + str(iQf)] = np.array( affine_components['f' + str(iQf)] )
                rhs_affine['f' + str(iQf)] = np.reshape( rhs_affine['f' + str(iQf)], \
                                                         rhs_affine['f' + str(iQf)].shape[0] )
            affine_components = rhs_affine
        
        return affine_components

    def assemble_fom_matrix( self, _param, _fom_specifics, _elements = [], _indices = []):
        
        if len( _elements ) == 0:
            matrix = self.M_engine.assemble_fom_matrix( self.convert_parameter( _param ), _fom_specifics )
            A = np.array( matrix['A'] )
            A[:, 0:2] = A[:, 0:2] - 1
            return A
        else:
            # if I'd convert elements and indices to int it also retrieve from int values inside the matrx from MATLAB
            # therefore I convert them to double
            matrix = self.M_engine.assemble_fom_matrix( self.convert_parameter( _param ), _fom_specifics, \
                                                        self.convert_parameter( _elements ), \
                                                        self.convert_parameter( _indices + 1 ) )
            
            A = np.array( matrix['A'] )
            A[:, 0:2] = A[:, 0:2] - 1
            
            return A
        
    # NB the +1 is needed to convert the python indices over MATLAB
    def find_mdeim_elements_fom_specifics( self, _fom_specifics, _indices_mat ):

        return np.array( self.M_engine.find_mdeim_elements_fom_specifics( _fom_specifics, \
                                       self.convert_indices( _indices_mat + 1 ) ) ).astype(int)
