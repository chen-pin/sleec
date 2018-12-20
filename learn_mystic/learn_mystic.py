"""Module for learning the Mystic optimization package"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import mystic as mys
from mystic import models 


class Quad(models.abstract_model.AbstractModel):
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
#        print("<evaluate(coeffs={}, coord={})".format(coeffs, coords))
        result = (coords[0]-coeffs[0])**2 + (coords[1]-coeffs[1])**2 
#        print("result = {}".format(result))
        return result

    def ForwardFactory(self, coeffs):
        """Requisite ForwardFactory method. 

        """
        def forward_quad(x):
            return self.evaluate(coeffs, x)
        return forward_quad


def _test():
    model = Quad()
    target = [0.0, 0.0]  # True values
    print("target = {}".format(target)) 
    num_data_pts = 1000
    coord_ampl = 5
    coords = coord_ampl*2.0*(np.random.rand(num_data_pts, 2)-0.5)
    print(coords.shape)
    p0 = 10*(np.random.rand(2)-0.5)  # Initial guess of coefficients
    print("Initial guess of coefficients = {}".format(p0))
    # Generate fake data
    noise_ampl = 0.1
    obsrvd_data = None
    obsrvd_data = np.array([(model.evaluate(coeffs=target, coords=coords[i,:]) 
                   + noise_ampl*2.0*(np.random.rand()-0.5)) for i in 
                   range(num_data_pts)])
    fig = plt.figure()
    plot = Axes3D(fig)
    plot.scatter(xs=coords[:, 0], ys=coords[:, 1], zs=obsrvd_data)
    print("starting point = {}".format(p0))
    #model.CostFactory = model.evaluate
    solution = mys.solvers.fmin_powell(
                   model.CostFactory2(coords.reshape(-1, num_data_pts)
                                      , obsrvd_data, 2), p0
                                      )
    print("solution = {}".format(solution))
    cost = model.__cost__(solution)
    print("solution cost value = {}".format(cost))
    predctn = np.array([model.evaluate(solution, coords=coords[i, :]) for i in 
               range(num_data_pts)])
    rms = np.sqrt(np.square((predctn-obsrvd_data)).mean())
    print("RMS predictn - data = {}".format(rms))
    plot.scatter(coords[:, 0], coords[:, 1], predctn)
    plt.show()

if __name__ == 'main':
    _test()
    
