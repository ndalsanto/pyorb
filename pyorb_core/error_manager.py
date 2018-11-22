#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 29 10:37:27 2018

@author: Niccolo" Dal Santo
@email : niccolo.dalsanto@epfl.ch
"""

error_switcher = {
        "Exception"  : Exception, \
        "SystemError": SystemError, \
        "ValueError" : ValueError, \
        "NameError"  : NameError, \
        "default"    : Warning
    }


def error_raiser( _error_type="Exception", _function="undefined_function", _message="This error has not been configured with any message, you should definitely do this" ):
    
    error_function = error_switcher.get( _error_type, SyntaxError )
    
    if error_function is SyntaxError:
        raise SyntaxError( "The chosen error '" + _error_type + "' is not available. The error message in any case is \n \n " + _message + " in function " + _function )
    
    raise error_function( _message + " in function " + _function )
    
    pass
