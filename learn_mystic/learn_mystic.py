"""Module for learning the Mystic optimization package"""

import numpy as np
import matplotlib.pyplot as plt
import mystic as mys


class Quad(mys.models.abstract_model.AbstractModel):
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
    def __init__(self, name='quad', metric=lambda x: np.sum(x*x)
                 , sigma=1.0 ):  
        super().__init__(name, metric, sigma)

    def evaluate(self, coeffs, coords):
        """Requisite model evaluation method.

        """
        print("<evaluate(coeffs={}, coord={})".format(coeffs, coords))
        result = ((coords[:, 0]-coeffs[0])**2 + (coords[:, 1]-coeffs[1])**2 
                + (coords[:, 2]-coeffs[2])**2)
        print("result = {}".format(result))
        return result

    def ForwardFactory(self, coeffs):
        """Requisite ForwardFactory method. 

        """
        def forward_quad(x):
            return self.evaluate(coeffs, x)
        return forward_quad


def _test():
    target = [0.0, 0.0, 0.0]  # True values
    print("target = {}".format(target)) 
    coords = np.random.rand(10, 3)
    p0 = 10*(np.random.rand(3)-0.5)  # Initial point
    print("starting point = {}".format(p0))
    model = Quad()
    #model.CostFactory = model.evaluate
    solution = mys.solvers.fmin_powell(model.CostFactory(target, coords), p0)
    print("solution = {}".format(solution))


if __name__ == 'main':
    _test()
    
