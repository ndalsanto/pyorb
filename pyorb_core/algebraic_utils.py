#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 15:42:29 2019

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

import numpy as np

# mat is a matrix in COO format, vec is a vector, such that the results is Av = mat * vec

def sparse_matrix_vector_mul( mat, vec ):

    Av = np.zeros( vec.shape  )

    nnz = mat.shape[0]

    for i in range( nnz ):
        Av[ int(mat[i, 0]), :] = Av[ int(mat[i, 0]), :] + mat[i, 2] * vec[ int(mat[i, 1]), : ]

    return Av
