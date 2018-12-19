"""Module for learning the Mystic optimization package"""

import numpy as np
import matplotlib.pyplot as plt
import mystic as mys


class Quad(mys.AbstractModel):
    """Toy Mystic model
    
    Parameters
    ----------
    a, b, c : float
        Coefficients

    Notes
    -----
    - Mystic module <mogi.py> is a good example of how to construct a mystic 
      base model.  
    
    """
    def __init__(self, name='quad', metric=lambda x: numpysum(x*x)
                 , sigma=1.0 ):  
        super().__init__(name, metric, sigma)

    def evaluate(self, coeffs, x):
        """Requisite model evaluation method.

        """
        return (x-coefs[0])**2 + (y-coefs[1])**2 + (z-coefs[2])**2

    def ForwardFactory(self, coeffs):
        """Requisite ForwardFactory method. 

        """
        def forward_quad(x):
            return self.evaluate(coeffs, x)
        return forward_quad


If __name__ == '__main__':
    
