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
import time

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

        return matlab.double( _np_array.astype(float).tolist() )

    def convert_indices( self, _indices ):

        return matlab.int64(_indices.tolist())

    def convert_types( self, _fom_specifics ):
        
        converted_fom_specifics = {}
        
        for key in _fom_specifics:
            
            if type(_fom_specifics[key]) == np.ndarray:
                converted_fom_specifics.update( { key : self.convert_double( _fom_specifics[key] ) } )
            else:
                converted_fom_specifics.update( { key : _fom_specifics[key] } )
            
        return converted_fom_specifics

    def solve_parameter( self, _param, _fom_specifics ):

        converted_fom_specifics = self.convert_types( _fom_specifics )

        _injected_param = self.convert_parameter( _param )
        
        u = self.M_engine.solve_parameter( _injected_param, converted_fom_specifics )

        sol = np.array( u['u'] )

        return sol[:, 0]

    # normally a MATLAB application can directly provide a dictionary with all the affine components
    def build_fom_affine_components( self, _operator, _num_affine_components, _fom_specifics ):

        print( 'Building affine components for operator %s' % _operator )
        
        converted_fom_specifics = self.convert_types( _fom_specifics )
        
        affine_components = self.M_engine.build_fom_affine_components( _operator, converted_fom_specifics )

        print( 'Finished to build affine components for operator %s ' % _operator )

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
        
        print('Returning affine components for operator %s' % _operator )

        return affine_components

    def assemble_fom_matrix( self, _param, _fom_specifics, _elements = [], _indices = []):
        
        _injected_param = self.convert_parameter( _param )

        converted_fom_specifics = self.convert_types( _fom_specifics )

        if len( _elements ) == 0:
            matrix = self.M_engine.assemble_fom_matrix( _injected_param, converted_fom_specifics )
            A = np.array( matrix['A'] )
            A[:, 0:2] = A[:, 0:2] - 1
            return A
        else:
            # if I'd convert elements and indices to int it also retrieve from int values inside the matrx from MATLAB
            # therefore I convert them to double
            matlab_elements = self.convert_parameter( _elements )
            matlab_indices  = self.convert_parameter( _indices + 1 )
            
            matrix = self.M_engine.assemble_fom_matrix( _injected_param, converted_fom_specifics, \
                                                        matlab_elements, matlab_indices )
            A = np.array( matrix['A'] )
            A[:, 0:2] = A[:, 0:2] - 1
            return A

    def assemble_fom_rhs( self, _param, _fom_specifics, _elements = [], _indices = []):

        converted_fom_specifics = self.convert_types( _fom_specifics )
        converted_params = self.convert_parameter( _param )
        converted_elements = self.convert_parameter( _elements )
        converted_indices = self.convert_parameter( _indices + 1 )

        if len( _elements ) == 0:
            rhs = self.M_engine.assemble_fom_rhs( self.convert_parameter( _param ), converted_fom_specifics )
            ff = np.array( rhs['f'] )
            ff = np.reshape( ff, (ff.shape[0], ) )
            return ff
        else:
            
            start = time.time()

            # if I convert elements and indices to int it would also retrieve from int values inside the matrx from MATLAB
            # therefore I convert them to double
            rhs = self.M_engine.assemble_fom_rhs( converted_params, converted_fom_specifics, \
                                                  converted_elements, \
                                                  converted_indices )

            end = time.time()
            print( 'Time to call rhs assembler from matlab' )
            print( end - start )

            ff = np.array( rhs['f'] )
#            self.M_engine.workspace["rhs_ptr"] = rhs['f']
#            ff = np.array( self.M_engine.eval("rhs_ptr.Value") )

            return ff
        
    # NB the +1 is needed to convert the python indices over MATLAB
    def find_deim_elements_fom_specifics( self, _fom_specifics, _indices ):

        converted_fom_specifics = self.convert_types( _fom_specifics )

        return np.array( self.M_engine.find_elements_for_deim_fom_specifics( converted_fom_specifics, \
                                       self.convert_indices( _indices + 1 ) ) ).astype(int)

    # NB the +1 is needed to convert the python indices over MATLAB
    def find_mdeim_elements_fom_specifics( self, _fom_specifics, _indices_mat ):

        converted_fom_specifics = self.convert_types( _fom_specifics )

        return np.array( self.M_engine.find_mdeim_elements_fom_specifics( converted_fom_specifics, \
                                       self.convert_indices( _indices_mat + 1 ) ) ).astype(int)


    def compute_natural_norm( self, _solution, _fom_specifics ):
        
        converted_fom_specifics = self.convert_types( _fom_specifics )

        return self.M_engine.compute_natural_norm( self.convert_double( _solution ), converted_fom_specifics )

    def assemble_fom_natural_norm_matrix( self, _fom_specifics ):
        
        converted_fom_specifics = self.convert_types( _fom_specifics )

        matrix = self.M_engine.assemble_fom_natural_norm_matrix( converted_fom_specifics )
        A = np.array( matrix['A'] )
        A[:, 0:2] = A[:, 0:2] - 1
        return A





