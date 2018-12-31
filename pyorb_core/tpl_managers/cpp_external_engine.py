#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 22 13:52:15 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import ctypes
from mpi4py import MPI
import numpy as np

from pyorb_core.tpl_managers import external_engine as ee

class c_fom_specifics( ctypes.Structure ):

    _fields_ = [('model', ctypes.POINTER( ctypes.c_char ) ), \
            ('datafile_path', ctypes.POINTER( ctypes.c_char ) ), \
            ('external_communicator', ctypes.c_void_p),\
            ('u', ctypes.POINTER( ctypes.c_double ) ),\
            ('A', ctypes.POINTER( ctypes.c_double ) ),\
            ('f', ctypes.POINTER( ctypes.c_double ) )]

class cpp_external_engine( ee.external_engine ):

    def __init__( self, _engine_type, _library_path ):

        ee.external_engine.__init__( self, _engine_type, _library_path )

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

        u = np.zeros( 0 ); cp_u = u.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
        A = np.zeros( 0 ); cp_A = A.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
        f = np.zeros( 0 ); cp_f = f.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )

        c_fom_spec = c_fom_specifics( ctypes.create_string_buffer( _fom_specifics['model'].encode('utf-8') ), \
                               ctypes.create_string_buffer( _fom_specifics['datafile_path'].encode('utf-8') ), \
                               ctypes.c_void_p( MPI._addressof( self.M_comm ) ), \
                               cp_u, cp_A, cp_f )

        compute_only_the_dim = True
        vector_dim = self.M_c_lib.solve_parameter( self.convert_parameter( _param ), c_fom_spec, compute_only_the_dim )

        u = np.zeros( vector_dim ); cp_u = u.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )

        c_fom_spec = c_fom_specifics( ctypes.create_string_buffer( _fom_specifics['model'].encode('utf-8') ), \
                               ctypes.create_string_buffer( _fom_specifics['datafile_path'].encode('utf-8') ), \
                               ctypes.c_void_p( MPI._addressof( self.M_comm ) ), \
                               cp_u, cp_A, cp_f )

        # the FOM library is supposed to fill in c_fom_specs.c_sol with the FOM
        compute_only_the_dim = False
        self.M_c_lib.solve_parameter( self.convert_parameter( _param ), c_fom_spec, compute_only_the_dim )

        return u

    def solve_parameter_and_export_solution( self, _param, _fom_specifics ):
        u = np.zeros( 0 ); cp_u = u.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
        A = np.zeros( 0 ); cp_A = A.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )
        f = np.zeros( 0 ); cp_f = f.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )

        c_fom_spec = c_fom_specifics( ctypes.create_string_buffer( _fom_specifics['model'].encode('utf-8') ), \
                               ctypes.create_string_buffer( _fom_specifics['datafile_path'].encode('utf-8') ), \
                               ctypes.c_void_p( MPI._addressof( self.M_comm ) ), \
                               cp_u, cp_A, cp_f )

        compute_only_the_dim = True
        vector_dim = self.M_c_lib.solve_parameter( self.convert_parameter( _param ), c_fom_spec, compute_only_the_dim )

        u = np.zeros( vector_dim ); cp_u = u.ctypes.data_as( ctypes.POINTER( ctypes.c_double ) )

        c_fom_spec = c_fom_specifics( ctypes.create_string_buffer( _fom_specifics['model'].encode('utf-8') ), \
                               ctypes.create_string_buffer( _fom_specifics['datafile_path'].encode('utf-8') ), \
                               ctypes.c_void_p( MPI._addressof( self.M_comm ) ), \
                               cp_u, cp_A, cp_f )

        self.M_c_lib.solve_parameter_and_export_solution( self.convert_parameter( _param ), c_fom_spec )

        return u

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

    M_comm = None
    M_c_lib = None
