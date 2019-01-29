"""Module for modeling a coronagraph using HCIPy as the propagation code."""

import os
import matplotlib.pyplot as plt
import pickle
import time
import numpy as np
import scipy.sparse as spr
import astropy.io.fits as fits
from hcipy import hcipy as hci


class TransmMask(hci.OpticalElement):
    """A `hci.OpticalElement` object representing a focal-plane mask.
    
    Parameters
    ----------
    input_grid : hci.Grid
        The grid on which the incoming wavefront is defined.
    mask_transm : ndarray
        Mask transmission function.  
    wavelength : float
        
    """
    def __init__(self, input_grid, mask_transm, wavelength=1):
        self.input_grid = input_grid
        self.mask_transm = mask_transm
        self.wavelength = wavelength

    def forward(self, wavefront):
        wf = wavefront.copy()
        wf.electric_field *= self.mask_transm
        return wf

    def backward(self, wavefront):
        wf = wavefront.copy()
        wf.electric_field /= self.mask_transm
        return wf


class FraunhoferBack(hci.OpticalElement):
    """`hcipy.OpticalElement` for backward Fraunhofer propagation.
    
    Parameters
    ----------
    pupil_grid : hci.Grid
        Pupil-plane grid
    focal_grid : hci.Grid
        Focal-plane grid
    wavelength_0 : float
        Reference wavelength
    focal_length : float
        Focal distance
    wavelength : float 
        Wavelength
    """

    def __init__(self, pupil_grid, focal_grid, wavelength_0=1, focal_length=1):
        self.pupil_grid = pupil_grid
        self.focal_grid = focal_grid
        self.wavelength_0 = wavelength_0
        self.focal_length = focal_length
        self.prop = hci.FraunhoferPropagator(input_grid=self.pupil_grid
                                             , output_grid=self.focal_grid
                                             , wavelength_0=self.wavelength_0
                                             , focal_length=self.focal_length
                                            )
    
    def forward(self, wavefront):
        return self.prop.backward(wavefront)
    def backward(self, wavefront):
        return self.prop.forward(wavefront)


class FresnelBack(hci.OpticalElement):
    """`hcipy.OpticalElement` for backward Fresnel propagation.

    Parameters
    ----------
    input_grid : hci.Grid
        Grid on which the input wavefront is defined.
    distance : float
        Propagation distance
    """
    def __init__(self, input_grid, distance):
        self.input_grid = input_grid
        self.distance = distance
        self.prop = hci.FresnelPropagator(input_grid=self.input_grid
                                         , distance=self.distance)
    def forward(self, wavefront):
        return self.prop.backward(wavefront)
    def backward(self, wavefront):
        return self.prop.forward(wavefront)


class InflncFtn(object): 
    """Make influence functions for `hci.DeformableMirror` using sparse 
    matrices.  
    
    Parameters
    ----------
    pupil_grid : hci.Grid
        The grid on which to calculate the influence functions spanning 
        the aperture of the DM.
    acttr_filnm : str
        Fits filename for a single actuator influence function (2darray).  
        The same function is applied to all actuators. 
    save_filnm : str
        Filename for storing/loading array of influence functions in NPZ 
        format.  
    num_acttr_across : int
        Number of actuators across the DM pupil.
    acttr_pxl_spacing: int
        Pixel scale (i.e. number of pixels between adjacent actuators) 
        associated with the single-actuator influence function (`self.acttr'). 
        For example, `self.acttr' loaded from file 
        'superinf-20180730-10x10.fits' is sampled at a scale of 10 pixels 
        between adjacent actuators.   
    
    Attributes
    ----------
    acttr_positns : hci.grid, @property
        Actuator center positions. 
    acttr_spacing : float, @property.
        Spacing between adjacent actuators.
    acttr : 2darray, @property
        Single-actuator influence function extracted from acttr_filnm
    acttr_grid : hci.grid, @property
        Automatically generated grid for a single actuator. Only used for 
        resampling/interpolation. 
    
    Notes
    -----
    - Future plan:  Expand to use a separate influence function for each 
      actuator.  
 
    """
    def __init__(self, pupil_grid=hci.make_pupil_grid(480)
                 , acttr_filnm=os.path.join(os.path.expanduser('~'), 'sleec_in'
                                            , 'superinf-20180730-10x10.fits')
                 , save_filnm=os.path.join(os.path.expanduser('~')
                                           , 'sleec_out', 'acttr3.npz')
                 , num_acttr_across=48, acttr_pxl_spacing=10):
        self.pupil_grid = pupil_grid
        self.acttr_filnm = acttr_filnm
        self.save_filnm = save_filnm 
        self.num_acttr_across = num_acttr_across
        self.acttr_pxl_spacing = acttr_pxl_spacing

    def generate(self, x_tilt=0.0, y_tilt=0.0, z_tilt=0.0, save=True):
        """Create DM influence functions.  If `save`=True, save the influence 
        functions as CSR  sparse matrices.  If not, return a hci.ModeBasis 
        object representing the array of influence functions. 

        Parameters
        ----------
        x_tilt : scalar
            Tilt of the DM around the x-axis in radians.
        y_tilt : scalar
            Tilt of the DM around the y-axis in radians.
        z_tilt : scalar
            Tilt of the DM around the z-axis in radians.
        save : boolean
            True --> save the array of influence functions as a CSR sparse 
            matrix in npz file format.  
            False --> return a hci.ModeBasis object representing the 
            array of influence functions.  

        Returns
        -------
        result : hci.ModeBasis, if `save`=True
            Influence functions 
        
        """
        evaluated_grid = (self.pupil_grid.scaled(1.0/np.cos([y_tilt, x_tilt]))
                         .rotated(-z_tilt))
        acttr_resample = hci.make_linear_interpolator_separated(
                         self.acttr.ravel(), self.acttr_grid, 0
                                                           )
        if save:
            modes = spr.lil_matrix((self.acttr_positns.size
                                   , evaluated_grid.size))
            for i in range(self.acttr_positns.size):
                modes[i, ] = acttr_resample(evaluated_grid
                                   .shifted(-self.acttr_positns[i]))
            spr.save_npz(self.save_filnm, modes.tocsr())
        else:
            modes =  [acttr_resample(evaluated_grid.shifted(-p)) for p in 
                      self.acttr_positns]
            modes = [hci.Field(m, self.pupil_grid) for m in modes]
            return hci.ModeBasis(modes)

    @property
    def acttr(self):
        """
        Returns
        -------
        acttr : 2darray, @property
            Single-actuator influence function extracted from acttr_filnm
        """
        return fits.getdata(self.acttr_filnm)

    @property
    def acttr_spacing(self):
        """
        Returns
        -------
        acttr_spacing : float
            Spacing between adjacent actuators.
        """
        return 1.0/self.num_acttr_across

    @property
    def acttr_positns(self):
        """
        Generate and automatic, uniform grid of actuator positions.  

        Returns
        -------
        A `hci.Grid` of actuator positions.
 
        """  
        extent = self.acttr_spacing*(self.num_acttr_across-1)
        return hci.make_uniform_grid(self.num_acttr_across, [extent]*2)

    @property
    def acttr_grid(self):
        """Generate an automatic grid for a single actuator.
        
        Notes
        -----
        - Only used for resampling/interpolation.  
        """
        return hci.make_uniform_grid(self.acttr.shape
                   , (self.acttr_spacing/self.acttr_pxl_spacing)
                     *np.array(self.acttr.shape)
                                    )
            
      
def _test(dm1_array=np.zeros(48*48), dm2_array=np.zeros(48*48)
                   , pupil_dia_pxl=256, fpm_wid=6.0, lyot_dia=0.9
                   , oversmpl=8, num_airy=24, dm_dia_num_acttr=48
                   , dm_seprtn=1/600e-6):
    """Test routine.
    
    """
    t0 = time.time()
    optical_elements = []
    # Make pupil
    pupil_grid = hci.make_pupil_grid(pupil_dia_pxl)
    focal_grid = hci.make_focal_grid(pupil_grid, oversmpl, num_airy)
    pupil_transm = hci.circular_aperture(1)(pupil_grid)
    pupil = TransmMask(pupil_grid, pupil_transm)
    optical_elements.append(pupil)
    # Influence Functions
# N.B. Re. Influence Fuctions:  
# An Hcipy influence function is an hcipy.ModeBasis.  A ModeBasis object is 
# a Python list comprising a list of Hcipy Field objects, each corresponding 
# to an actuator.   
    pickle_fn = ('inflnc_ftns_'+str(dm_dia_num_acttr)+'_'+str(pupil_dia_pxl)
                 +'.pkl')
    pickle_path = os.path.join(os.path.expanduser('~'), 'data_sleec', pickle_fn)
    try:
        f = open(pickle_path, 'rb')
    except FileNotFoundError:
        print("Influence-functions pickle file not found.")
        print("Generate influence functions on the fly.")
        inflnc_ftns = hci.make_xinetics_influence_functions(pupil_grid
                                                       , dm_dia_num_acttr
                                                       , 1.0/dm_dia_num_acttr)
        f = open(pickle_path, 'wb')
        pickle.dump(inflnc_ftns, f)
    else:
        inflnc_ftns = pickle.load(f)
        inflnc_ftns[0] = hci.Field(inflnc_ftns[0], pupil_grid)
    # DM1
    dm1 = hci.DeformableMirror(inflnc_ftns)
#    dm1.actuators = np.zeros(len(inflnc_ftns))
#    dm1.actuators = np.sqrt(1e-3)*np.random.rand(len(inflnc_ftns))
    dm1.actuators = dm1_array
    optical_elements.append(dm1)
    # Propagator to DM2
    prop_to_dm2 = hci.FresnelPropagator(pupil_grid, dm_seprtn)
    optical_elements.append(prop_to_dm2)
    # DM2
    dm2 = hci.DeformableMirror(inflnc_ftns)
#    dm2.actuators = np.zeros(len(inflnc_ftns))
#    dm2.actuators = np.sqrt(1e-3)*np.random.rand(len(inflnc_ftns))
    dm2.actuators = dm2_array
    optical_elements.append(dm2)
    # Propagator to Pupil
    prop_to_pupil = FresnelBack(pupil_grid, dm_seprtn)
    optical_elements.append(prop_to_pupil)
    # Propagator to Focal Plane
    prop_to_focal = hci.FraunhoferPropagator(pupil_grid, focal_grid)
    optical_elements.append(prop_to_focal)
    # Focal-Plane Mask
    fpm_opacity = hci.circular_aperture(fpm_wid)
    fpm_transm = 1.0 - fpm_opacity(focal_grid)
    fpm = TransmMask(focal_grid, fpm_transm)
    optical_elements.append(fpm)
    # Propagator to Lyot Plane
    prop_to_lyot = FraunhoferBack(pupil_grid, focal_grid)
    optical_elements.append(prop_to_lyot)
    # Lyot-Plane Mask
    lyot_transm = hci.circular_aperture(lyot_dia)(pupil_grid)
    lyot_mask = TransmMask(pupil_grid, lyot_transm)
    optical_elements.append(lyot_mask)
    # Propagator to Camera
    prop_to_camera = hci.FraunhoferPropagator(pupil_grid, focal_grid)
    optical_elements.append(prop_to_camera)
    # Generate Non-Corongraphic PSF
    psf = prop_to_focal(hci.Wavefront(pupil_transm))
    coro = hci.OpticalSystem(optical_elements)
    camera_wf = coro(hci.Wavefront(pupil_transm))
    print("Computation took {:g} s".format(time.time()-t0))
    hci.imshow_field(np.log10(camera_wf.intensity/psf.intensity.max())
                              , vmin=-12, cmap='gray')
    plt.colorbar()
    plt.show()
#    plt.savefig('dm2_int.jpg')
    plt.close()
    hci.imshow_field(camera_wf.phase, vmin=-1, vmax=1, cmap='gray')
#    plt.colorbar()
    plt.show()
#    plt.savefig('dm2_ph.jpg')
    plt.close()
#    reshape_wid = int(np.sqrt(camera_wf.phase.shape[0]))
#    diff = (camera_wf.phase.reshape((reshape_wid, -1))[reshape_wid//2, :]
#           - camera_wf.phase.reshape((reshape_wid, -1))[:, reshape_wid//2])
#    plt.plot(diff)
#    plt.show()


if __name__ == '__main__':
    _test()
