"""
Module for generating Slepian functions using Frederik Simons' code. 

Notes
-----
This module depends on the following packages: 
- Python packages 
    - numpy
    - oct2py
    - hcipy
- GNU software:  Octave.  `Ocatve` needs to be installed on the system.   
- Frederik Simons' Matlab scripts, downloadable from `Simons Slepian 
<https://github.com/harig00/MMASS>'_.  Only the folowing scripts are required. 
Either place the `.m` files in the working directory (where python is run), or 
use oct2py.octave.addpath('path') to add the path where the `.m` files reside. 
`swdisk.m` is the main/frontend script.   
    - defval.m
    - devilliers.m
    - difer.m
    - gausslegendrecof.m
    - ht.m
    - indeks.m
    - jacobi.m
    - kernelc2D.m
    - kernelc2D_pin.m
    - legendrecof.m
    - legendreroots.m
    - localization2D.m
    - nanmean.m
    - nsubm.m
    - rindeks.m
    - star69.m
    - swdisk.m
    - swfried2d.m
    - swspace2d.m
    - swspectral2d.m
    - swvals2d.m
    - tridiag.m

Reference
---------
simons, F.J., Wang, D.V. (2011), "Spatiospectral concentration in the 
Cartesian plane," *Int J Geomath*, DOI 10.1007/s13137-011-0016-z."  

"""
import numpy as np
from oct2py import octave as oct
import hcipy as hci


def gen_slepians_r(n2d, L, method, orthgnlty, D, grid, num_dense_pts=4096):
    """
    Return a mode basis of radial Slepian functions concentrated to a 
    circular domain.  Assume m = 0 and k=1.  Slepians corresponding to all 
    elements of 'n2d' is computed.  See documentation in <swdisk.m> by 
    F. Simons for details.  

    Parameters
    ----------
    n2d : array-like or integer
        The two-dimensional Shannon number.  For a Fraunhofer propagation, 
        n2d = (pi*N/4)^2, where N is the number of cycles across the 
        aperture. For example, if we are solving for a pupil-apodizing 
        function to limit energy in a 3wvl/D (radius) zone, n2d = 5.6.   
    L : The maximum degree in the Jacobi/Bessel expansion (default:  2*n2d),  
        or else the maximum polynomial degree for which GL wil be exact. 
        Simons states that increasing `L` is almost always a good idea,
        presumably because it increases accuracy.  However, P. Chen found that 
        computations can fail for high `L` values.
    method : str
        'DV' compute by the method of De Villiers
        'SE' by Slepian "extension" [the default method]
        'GL' by direct Gauss-Legendre integration
    orthgnlty : int
        1 Scales the solution for weightless orthogonality
        0 Scales the solution for orthogonality with weight r
    D : float
        Optical aperture diameter
    grid : hcipy.Grid
        The hcipy `Grid` object supporting the output function values. 
    num_dense_pts : int
        Dense sampling points for computing the Slepians.  Higher order 
        functions require dense sampling due to high order Bessel functions 
        used in the computation.   
    
    Returns
    -------
    return : `hci.ModeBasis`
        A list of the user specified Slepian functions evaluated on the input 
        grid (interpolated from the densely sampled computation). 

    Notes
    -----
    An `hcipy.ModeBasis` object is a list of `hcipy.Field` objects.  

    """
    polar_grid = grid.as_('polar')
    r = 2.0*polar_grid.r/D 
    # Compute the Slepian functions on a dense grid, which will then be used 
    # for interpolation onto the input grid.   
    r_dense_min = 0.0
    r_dense_max = 2.0 
    r_dense = np.linspace(r_dense_min, r_dense_max, num_dense_pts)
    modes = []
    for n2d_i in n2d:
        slp_dense = oct.swdisk(0, n2d_i, 1, L, r_dense, method, orthgnlty)
        modes.append(hci.Field(np.interp(r, r_dense, slp_dense[:, 0]), grid))
    return hci.ModeBasis(modes, grid)


def gen_slepians_r_k(m, n2d, k, L, method, orthgnlty, D, grid
                   , num_dense_pts=4096):
    """
    Return a mode basis of `k` radial Slepian functions for fixed, 
    user-specified `m` and `n2d` values.  See 
    documentation in <swdisk.m> by F. Simons for details.  

    Parameters
    ----------
    m : int
        The fixed order of the eigenvalue problem.  It is analogous to the 
        subscript n for Bessel functions J_n.
    n2d : float
        The two-dimensional Shannon number.  For a Fraunhofer propagation, 
        n2d = (pi*N/4)^2, where N is the number of cycles across the 
        aperture. For example, if we are solving for a pupil-apodizing 
        function to limit energy in a 3wvl/D (radius) zone, n2d = 5.6.   
    k : int
        The number of functions to output.  For example, k=3, the output 
        contains 3 Slepain functions g_1, g_2, g_3.  The g_k functions have 
        increasing frequencies with increasing k.  Examining a small number 
        examples suggest the number of zeros in the range [0, 1] equals k-1. 
        Note that k > n2d would not make much sense because the Shannon number  
        signifies the number of significant eigenfunctions (Slepians) that 
        can simultaneously be concentrated into the defined spatio-spectral 
        domain.
    L : The maximum degree in the Jacobi/Bessel expansion (default:  2*n2d),  
        or else the maximum polynomial degree for which GL wil be exact. 
        Simons states that increasing `L` is almost always a good idea,
        presumably because it increases accuracy.  However, P. Chen found that 
        computations can fail for high `L` values.
    method : str
        'DV' compute by the method of De Villiers
        'SE' by Slepian "extension" [the default method]
        'GL' by direct Gauss-Legendre integration
    orthgnlty : int
        1 Scales the solution for weightless orthogonality
        0 Scales the solution for orthogonality with weight r
    D : float
        Optical aperture diameter
    grid : hcipy.Grid
        The hcipy `Grid` object supporting the output function values. 
    num_dense_pts : int
        Dense sampling points for computing the Slepians.  Higher order 
        functions require dense sampling due to high order Bessel functions 
        used in the computation.   
    
    Returns
    -------
    return : `hci.ModeBasis`
        A list of the user specified Slepian functions evaluated on the input 
        grid (interpolated from the densely sampled computation). 
    Notes
    -----
    An `hcipy.ModeBasis` object is a list of `hcipy.Field` objects.  
    """
    polar_grid = grid.as_('polar')
    r = 2.0*polar_grid.r/D 
    # Compute the Slepian functions on a dense grid, which will then be used 
    # for interpolation onto the input grid.   
    r_dense_min = 0.0
    r_dense_max = 2.0 
    r_dense = np.linspace(r_dense_min, r_dense_max, num_dense_pts)
    slp_dense = oct.swdisk(m, n2d, k, L, r_dense, method, orthgnlty)
    modes = [hci.Field(np.interp(r, r_dense, slp_dense[:, i]), grid) 
                   for i in range(k)]
    return hci.ModeBasis(modes, grid)


def gen_slepian_fourier(m_max, n2d, k, L, method, orthgnlty, D, grid
        , num_dense_pts=4096): 
    """
    Return a mode basis of `k` Slepian-Fourier (radial-azimuthal) functions 
    for fixed, user-specified `m` and `n2d` values.  See 
    documentation in <swdisk.m> by F. Simons for details.  

    Parameters
    ----------
    m_max : int
        Maximum absolute `m` value.  A list of `m` values from -`m_max` to 
        `m_max` will be generated.  `m` is the fixed order of the eigenvalue 
        problem.  It is analogous to the subscript n for Bessel functions J_n.
    n2d : float
        The two-dimensional Shannon number.  For a Fraunhofer propagation, 
        n2d = (pi*N/4)^2, where N is the number of cycles across the 
        aperture. For example, if we are solving for a pupil-apodizing 
        function to limit energy in a 3wvl/D (radius) zone, n2d = 5.6.   
    k : int
        The number of functions to output.  For example, k=3, the output 
        contains 3 Slepain functions g_1, g_2, g_3.  The g_k functions have 
        increasing frequencies with increasing k.  Examining a small number 
        examples suggest the number of zeros in the range [0, 1] equals k-1. 
        Note that k > n2d would not make much sense because the Shannon number  
        signifies the number of significant eigenfunctions (Slepians) that 
        can simultaneously be concentrated into the defined spatio-spectral 
        domain.
    L : The maximum degree in the Jacobi/Bessel expansion (default:  2*n2d),  
        or else the maximum polynomial degree for which GL wil be exact. 
        Simons states that increasing `L` is almost always a good idea,
        presumably because it increases accuracy.  However, P. Chen found that 
        computations can fail for high `L` values.
    method : str
        'DV' compute by the method of De Villiers
        'SE' by Slepian "extension" [the default method]
        'GL' by direct Gauss-Legendre integration
    orthgnlty : int
        1 Scales the solution for weightless orthogonality
        0 Scales the solution for orthogonality with weight r
    D : float
        Optical aperture diameter
    grid : hcipy.Grid
        The hcipy `Grid` object supporting the output function values. 
    num_dense_pts : int
        Dense sampling points for computing the Slepians.  Higher order 
        functions require dense sampling due to high order Bessel functions 
        used in the computation.   
    
    Returns
    -------
    return : `hci.ModeBasis`
        A list of the user specified Slepian-Fourier  functions evaluated on 
        the input grid (interpolated from the densely sampled computation). 
    Notes
    -----
    An `hcipy.ModeBasis` object is a list of `hcipy.Field` objects.  
    """
    polar_grid = grid.as_('polar')
    r = 2.0*polar_grid.r/D 
    theta = polar_grid.theta
    # Compute the Slepian functions on a dense grid, which will then be used 
    # for interpolation onto the input grid.   
    r_dense_min = 0.0
    r_dense_max = 2.0 
    r_dense = np.linspace(r_dense_min, r_dense_max, num_dense_pts)
    modes = []
    m = range(-m_max, m_max+1)
    for m_i in m:
        if m_i < 0:
            fourier = np.sin(np.abs(m_i)*theta)
        elif m_i == 0:
            fourier = 1.0
        else:
            fourier = np.cos(m_i*theta)
        for n2d_i in n2d:
            print(m_i, n2d_i)
            slp_dense = oct.swdisk(int(np.abs(m_i)), n2d_i, k, L, r_dense
                                   , method, orthgnlty)
            for i in range(k):
                modes.append(
                    hci.Field(np.interp(r, r_dense, slp_dense[:, i])*fourier
                                   , grid)
                             ) 
#            modes.append(
#                    [hci.Field(np.interp(r, r_dense, slp_dense[:, i])*fourier
#                               , grid) 
#                         for i in range(k)]
#                        )
    return hci.ModeBasis(modes, grid)
    

def _test_gen_slepian_fourier():
    import matplotlib.pyplot as plt
    m_max = 2
    n2d = [10, 20]
    k = 2
    L = 40
    method = 'SE'
    orthgnlty = 0
    D = 8.0
    grid = hci.make_pupil_grid(256, D)
    basis = gen_slepian_fourier(m_max, n2d, k, L, method, orthgnlty, D, grid
                   , num_dense_pts=4096)
    print("length of basis = {}".format(len(basis)))
    print("length of first element of basis = {}".format(len(basis[0])))
    print("first element of basis:  {}".format(basis[0]))
    for mode in basis:
        hci.imshow_field(mode)
        plt.colorbar()
        plt.show()
        plt.close()


if __name__ == '__main__':
    _test_gen_slepian_fourier()




