"""Module for learning the Mystic optimization package"""

import numpy as np
import matplotlib.pyplot as plt
import mystic as mys


def quad_ftn(x, y, z):
    """Toy function for learning mystic.
    
    Parameters
    ----------
    x, y, z : float
        Coordinate values
    
    Returns
    -------
    3-D quadratic function with minimum at the origin.

    """
    
    return x**2 + y**2 + z**2
    

class Quad(mys.AbstractModel):
    """Toy Mystic model
    
    Parameters
    ----------
    a, b, c : float
        Coefficients
    
    """

    def __init__(self, name='quad', metric=lambda x: numpysum(x*x)
                 , sigma=1.0 ):  
        super().__init__(name, metric, sigma)


    def evaluate(self, x, y, z):
        """Requisite model evaluation method.

        """

        return quad_ftn(x, y, z)

    def ForwardFactory(self, x, y, z):
        """Requisite ForwardFactory method. 

        """

        return self.evaluate


If __name__ == '__main__':
    
