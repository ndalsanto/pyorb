#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 10 10:24:47 2018

@author: Niccolo' Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

#%%

import numpy as np

import sys
sys.path.insert(0, '../../core')

import rb_manager as rm

print( rm.__doc__ )


my_rb_manager = rm.RbManager( )

#my_rb_manager.set_random_snapshots_matrix( 4, 2 )
#my_rb_manager.print_snapshots_matrix( )
#my_rb_manager.perform_pod( )






snapshots_file = 'train_snapshots_matrix_20_50.txt'

my_rb_manager.import_snapshots_matrix( snapshots_file )
print( "The new number of snapshots is %d " % my_rb_manager.get_number_of_snapshots( ) )
my_rb_manager.perform_pod( )



