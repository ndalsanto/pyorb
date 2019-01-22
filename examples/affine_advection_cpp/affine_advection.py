# @Author: Luca Pegolotti
# @Date:   2018-12-18T16:59:45+01:00
# @Email:  luca.pegolotti@epfl.ch
# @Last modified by:   Luca Pegolotti
# @Last modified time: 2018-12-18T17:45:47+01:00


import numpy as np
import math

import pyorb_core.pde_problem.fom_problem as fp

def aa_theta_a( _param, _q ):
    assert _q <= 2

    theta = _param[1] * math.pi / 180.0
    if( _q == 0 ):
        return _param[0]
    elif( _q == 1 ):
        return math.sin( theta )
    elif( _q == 2 ):
        return math.cos( theta )


def aa_theta_f( _param, _q ):
    return 1.0


class affine_advection_problem( fp.fom_problem ):

    def __init__( self, _parameter_handler ):
        fp.fom_problem.__init__( self, _parameter_handler )

        return

    def define_theta_functions( self ):

        self.M_theta_a = aa_theta_a
        self.M_theta_f = aa_theta_f

        return
