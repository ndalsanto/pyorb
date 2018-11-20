#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 25 14:17:29 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import matlab.engine
import ctypes
import numpy as np
from mpi4py import MPI

import error_manager as em


class external_engine( ):

    def __init__( self, _engine_type, _library_path ):
        
        self.M_engine_type  = _engine_type
        self.M_library_path = _library_path
        
        return

    def start_engine( self ):
        self.start_specific_engine( )
        return
    
    def quit_engine( self ):
        self.quit_specific_engine( )
        return

    def start_specific_engine( self ):
        
        em.error_raiser( 'SystemError', 'external_engine::start_engine', "You are using the default start_specific_engine, \
                          please provide specific ones for your specific engine " )
        return

    def quit_specific_engine( self ):
        
        em.error_raiser( 'SystemError', 'external_engine::quit_engine', "You are using the default quit_specific_engine, \
                          please provide specific ones for your specific engine " )
        return

    def convert_parameter( self, _param ):

        em.error_raiser( 'SystemError', 'external_engine::convert_parameter', "You are using the default convert_parameter, \
                          please provide specific ones for your specific engine " )
        return

    def solve_parameter( self, _param, _fom_specifics ):

        em.error_raiser( 'SystemError', 'external_engine::solve_parameter', "You are using the default solve_parameter, \
                          please provide specific ones for your specific engine " )
        return

    def build_rb_affine_component( self, _param, _fom_specifics ):

        em.error_raiser( 'SystemError', 'external_engine::build_rb_affine_component', "You are using the default build_rb_affine_component, \
                          please provide specific ones for your specific engine " )
        return

    def build_fom_affine_components( self, _operator, _num_affine_components, _fom_specifics ):
        
        em.error_raiser( 'SystemError', 'external_engine::build_rb_affine_component', "You are using the default build_fom_affine_components, \
                          please provide specific ones for your specific engine " )
        return

    def assemble_fom_matrix( self, _param, _fom_specifics, _elements, _indices ):
        
        em.error_raiser( 'SystemError', 'external_engine::assemble_fom_matrix', "You are using the default assemble_fom_matrix, \
                          please provide specific ones for your specific engine " )
        return

    def find_mdeim_elements_fom_specifics( self, _fom_specifics, _indices_mat ):

        em.error_raiser( 'SystemError', 'external_engine::find_mdeim_elements_fom_specifics', \
                         "You are using the default find_mdeim_elements_fom_specifics, \
                          please provide specific ones for your specific engine " )
        return


    M_engine_type = ""
    M_library_path = ""
    M_engine = 0




class c_fom_specifics( ctypes.Structure ):

    _fields_ = [('model', ctypes.POINTER( ctypes.c_char ) ), \
            ('datafile_path', ctypes.POINTER( ctypes.c_char ) ), \
            ('external_communicator', ctypes.c_void_p),\
            ('u', ctypes.POINTER( ctypes.c_double ) ),\
            ('A', ctypes.POINTER( ctypes.c_double ) ),\
            ('f', ctypes.POINTER( ctypes.c_double ) )]




class cpp_external_engine( external_engine ):
    
    def __init__( self, _engine_type, _library_path ):
        
        external_engine.__init__( self, _engine_type, _library_path )
        
        return

    def start_specific_engine( self ):

        self.M_comm = MPI.COMM_WORLD
        self.M_c_lib = ctypes.cdll.LoadLibrary( self.M_library_path )
        
        return

    def quit_specific_engine( self ):

        return

    def convert_parameter( self, _param ):

        return _param.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )

    def convert_double( self, _np_array ):

        return _np_array.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )

    def convert_indices( self, _indices ):

        return _indices.ctypes.data_as( ctypes.POINTER( ctypes.c_int ) )

    def solve_parameter( self, _param, _fom_specifics ):
        
        u = np.zeros( 6267 )
        cp_u = u.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
        A = np.zeros( 0 )
        cp_A = A.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
        f = np.zeros( 0 )
        cp_f = f.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
        
        c_fom_spec = c_fom_specifics( ctypes.create_string_buffer( _fom_specifics['model'].encode('utf-8') ), \
                               ctypes.create_string_buffer( _fom_specifics['datafile_path'].encode('utf-8') ), \
                               ctypes.c_void_p( MPI._addressof( self.M_comm ) ), \
                               cp_u, cp_A, cp_f )

        # the FOM is supposed to fill in c_fom_specs.c_sol with the FOM 
        self.M_c_lib.solve_parameter( self.convert_parameter( _param ), c_fom_spec )

        sol = {}
        sol['u'] = u

        return sol

#    def build_rb_affine_component( self, _basis, _q, _operator, _fom_specifics ):
#
#        return self.M_engine.build_rb_affine_component( _basis, _q, _operator, _fom_specifics )


    def build_fom_affine_components( self, _operator, _num_affine_components, _fom_specifics ):
        
        if _operator == 'A':
            return self.build_fom_matrix_affine_components( _num_affine_components, _fom_specifics )
        elif _operator == 'f':
            return self.build_fom_vector_affine_components( _num_affine_components, _fom_specifics )
        else:
            em.error_raiser( 'SystemError', 'cpp_external_engine::build_fom_affine_components', \
                             "You should provide a valid operator for building the corresponding affine component" )
            

    
    def build_fom_matrix_affine_components( self, _num_affine_components, _fom_specifics ):

        # fake parameters which are not used for real
        u = np.zeros( 0 )
        cp_u = u.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
        f = np.zeros( 0 )
        cp_f = f.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
        
        affine_components = { }
        
        print( 'Building %d matrix affine components' % _num_affine_components )

        for iQa in range( _num_affine_components ):
            
            print( 'Building rhs affine components number %d ' % iQa )

            A = np.zeros( 0 )
            cp_A = A.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
            
            c_fom_spec = c_fom_specifics( ctypes.create_string_buffer( _fom_specifics['model'].encode('utf-8') ), \
                                          ctypes.create_string_buffer( _fom_specifics['datafile_path'].encode('utf-8') ), \
                                          ctypes.c_void_p( MPI._addressof( self.M_comm ) ), \
                                          cp_u, cp_A, cp_f )

            # retrieve at first the number of nnz
            compute_only_the_nnz = True
            matrix_nnz = self.M_c_lib.build_fom_affine_components( ctypes.create_string_buffer( 'A'.encode('utf-8') ), iQa, 
                                                                c_fom_spec, compute_only_the_nnz )
            
            # For the COO format, 3 * nnz must be stored
            A = np.zeros( matrix_nnz * 3 )
            cp_A = A.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
            c_fom_spec = c_fom_specifics( ctypes.create_string_buffer( _fom_specifics['model'].encode('utf-8') ), \
                              ctypes.create_string_buffer( _fom_specifics['datafile_path'].encode('utf-8') ), \
                              ctypes.c_void_p( MPI._addressof( self.M_comm ) ), \
                              cp_u, cp_A, cp_f )

            # retrieve at first the number of nnz
            compute_only_the_nnz = False
            self.M_c_lib.build_fom_affine_components( ctypes.create_string_buffer( 'A'.encode('utf-8') ), iQa, 
                                                      c_fom_spec, compute_only_the_nnz )

            affine_components['A' + str(iQa)] = np.reshape( A, (matrix_nnz, 3), order='C' )
        
        return affine_components


    def build_fom_vector_affine_components( self, _num_affine_components, _fom_specifics ):

        # fake parameters which are not used in the computation
        u = np.zeros( 0 )
        cp_u = u.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
        A = np.zeros( 0 )
        cp_A = A.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
        
        affine_components = { }
        
        print( 'Building %d rhs affine components' % _num_affine_components )
        
        for iQf in range( _num_affine_components ):

            print( 'Building rhs affine components number %d ' % iQf )

            f = np.zeros( 0 )
            cp_f = f.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
            
            c_fom_spec = c_fom_specifics( ctypes.create_string_buffer( _fom_specifics['model'].encode('utf-8') ), \
                                          ctypes.create_string_buffer( _fom_specifics['datafile_path'].encode('utf-8') ), \
                                          ctypes.c_void_p( MPI._addressof( self.M_comm ) ), \
                                          cp_u, cp_A, cp_f )

            # retrieve at first the dim of the rhs
            compute_only_the_dim = True
            vector_dim = self.M_c_lib.build_fom_affine_components( ctypes.create_string_buffer( 'f'.encode('utf-8') ), iQf, 
                                                                  c_fom_spec, compute_only_the_dim )
            
            f = np.zeros( vector_dim )
            cp_f = f.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
            
            c_fom_spec = c_fom_specifics( ctypes.create_string_buffer( _fom_specifics['model'].encode('utf-8') ), \
                                          ctypes.create_string_buffer( _fom_specifics['datafile_path'].encode('utf-8') ), \
                                          ctypes.c_void_p( MPI._addressof( self.M_comm ) ), \
                                          cp_u, cp_A, cp_f )
            
            compute_only_the_dim = False
            self.M_c_lib.build_fom_affine_components( ctypes.create_string_buffer( 'f'.encode('utf-8') ), iQf, 
                                                     c_fom_spec, compute_only_the_dim )
            
            affine_components['f' + str(iQf)] = f

        return affine_components


#    def assemble_fom_matrix( self, _param, _fom_specifics, _elements = [], _indices = []):
#        
#        if len( _elements ) == 0:
#            return self.M_engine.assemble_fom_matrix( self.convert_parameter( _param ), _fom_specifics )
#        else:
#            matrix = self.M_engine.assemble_fom_matrix( self.convert_parameter( _param ), _fom_specifics, \
#                                                        self.convert_parameter(_elements), \
#                                                        self.convert_parameter(_indices) )
#            return np.array( matrix['A'] )
#        
#    def find_mdeim_elements_fom_specifics( self, _fom_specifics, _indices_mat ):
#
#        return np.array( self.M_engine.find_mdeim_elements_fom_specifics( _fom_specifics, \
#                                       self.convert_indices( _indices_mat ) ) ).astype(int)
        

    M_comm = 0
    M_c_lib = 0




class matlab_external_engine( external_engine ):
    
    def __init__( self, _engine_type, _library_path ):
        
        external_engine.__init__( self, _engine_type, _library_path )
        
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

        return self.M_engine.solve_parameter( self.convert_parameter( _param ), _fom_specifics )

    def build_rb_affine_component( self, _basis, _q, _operator, _fom_specifics ):

        return self.M_engine.build_rb_affine_component( _basis, _q, _operator, _fom_specifics )

    # normally a MATLAB application can directly provide a dictionary with all the affine components 
    def build_fom_affine_components( self, _operator, _num_affine_components, _fom_specifics ):
        
        affine_components = self.M_engine.build_fom_affine_components( _operator, _fom_specifics )

        # rescale the matrix indices so that the counting starts from 0 (and not from 1 as in MATLAB)
        if _operator == 'A':
            for iQa in range( _num_affine_components ):
                affine_components['A' + str(iQa)][:, 0:2] = affine_components['A' + str(iQa)][:, 0:2] - 1

        return affine_components

    def assemble_fom_matrix( self, _param, _fom_specifics, _elements = [], _indices = []):
        
        if len( _elements ) == 0:
            return self.M_engine.assemble_fom_matrix( self.convert_parameter( _param ), _fom_specifics )
        else:
            # if I'd convert elements and indices to int it also retrieve from int values inside the matrx from MATLAB
            # therefore I convert them to double
            matrix = self.M_engine.assemble_fom_matrix( self.convert_parameter( _param ), _fom_specifics, \
                                                        self.convert_parameter(_elements), \
                                                        self.convert_parameter(_indices) )
            return np.array( matrix['A'] )
        
    def find_mdeim_elements_fom_specifics( self, _fom_specifics, _indices_mat ):

        return np.array( self.M_engine.find_mdeim_elements_fom_specifics( _fom_specifics, \
                                       self.convert_indices( _indices_mat ) ) ).astype(int)




class external_engine_manager( ):
    
    def __init__( self, _engine_type, _library_path ):
        
        self.M_engine_type  = _engine_type
        self.M_library_path = _library_path
        
        if _engine_type == 'matlab':
            self.M_external_engine = matlab_external_engine( _engine_type, _library_path )
        elif _engine_type == 'cpp':
            self.M_external_engine = cpp_external_engine( _engine_type, _library_path )
        
        return

    M_engine_type = ""
    M_library_path = ""
    M_external_engine = 0
    
    def get_external_engine( self ):
        return self.M_external_engine

    def start_engine( self ):
        self.M_external_engine.start_engine( )
        return
    
    def quit_engine( self ):
        self.M_external_engine.quit_engine( )
        return



