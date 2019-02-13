#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 22 17:02:27 2019

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""


def build_matlab_default_fom_specifics( ):
    
    fom_specifics =  { \
                'polynomial_degree' : 'P1', \
                'model'             : 'undefined_model', \
                'use_nonhomogeneous_dirichlet' : 'N', \
                'mesh_name' : '' }
    
    return fom_specifics
