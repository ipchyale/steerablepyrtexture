#!/usr/bin/python
"""(functions that interact with the C code, for handling convolutions mostly.)

REPLACED:
Python implementations of the C code functions. In this case only one is needed
for steerable pyramids in frequency space. (was pointOp in c)


"""
import numpy as np


def point_op(im, lut, origin, increment, warnings=False):
    """Apply a point operation, specified by lookup table `lut`, to `image`

    Arguments
    ---------
    image : `array_like`
        1d or 2d array
    lut : `array_like`
        a row or column vector, assumed to contain (equi-spaced) samples of the function.
    origin : `float`
        specifies the abscissa associated with the first sample
    increment : `float`
        specifies the spacing between samples.
    warnings : `bool`
        whether to print a warning whenever the lookup table is extrapolated

    """
    result = np.empty_like(im)
    # this way we can use python booleans when calling
    if warnings:
        warnings = 1
    else:
        warnings = 0
    X = origin + increment*np.arange(0,np.size(lut))
    Y = lut
    
    sz1 = np.size(im,0)
    sz2 = np.size(im,1)

    result = np.interp(im.flatten(), X, Y)

    result = result.reshape(sz1,sz2)

    return np.array(result)
