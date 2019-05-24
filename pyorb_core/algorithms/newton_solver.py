# -*- coding: utf-8 -*-

import numpy as np

def newton_solver( _jacobian, _residual, _u0, _tol, _n_max ):
    
    it = 0
    u_n = _u0
    
    res = _residual( u_n )
    res_norm = np.linalg.norm( res )
#    print( 'Starting residual norm is %e' % res_norm )
    
    while res_norm > _tol and it < _n_max:
        it = it + 1
        u_n_1 = 1.* u_n
        u_n = u_n_1 - np.linalg.solve( _jacobian(u_n_1), _residual(u_n_1) )
        res = _residual( u_n )

#        print( 'Residual is' )
#        print( res )

        res_norm = np.linalg.norm( res )
#        print( 'Iteration %d residual norm is %e' % (it, res_norm) )
    
    return u_n