"""Module for modeling a coronagraph using HCIPy as the propagation code."""

import os
import pickle
import time
import numpy as np
import scipy.sparse as spr
import scipy.interpolate as scint
import astropy.io.fits as fits
from hcipy import hcipy as hci


class DefrmblMirror(hci.OpticalElement):
    """A deformable mirror (DM) with input surface. This is a modification of 
    hci.DeformableMirror.  

    Parameters
    ----------
    surfc : hci.Field
        The DM surface.

    Attributes
    ----------
    phase_for : hci.Field, @ property
        DM surface in terms of phase. 

    """
    def __init__(self, surfc):
        self.surfc = surfc

    def forward(self, wavefront):
        """Propagate a wavefront through the deformable mirror.

        Parameters
        ----------
        wavefront : Wavefront
        The incoming wavefront.
        
        Returns
        -------
        Wavefront
        The reflected wavefront.
        """
        wf = wavefront.copy()
        wf.electric_field *= np.exp(
                                 2j*self.surfc*wavefront.wavenumber
                                   )
        return wf

    def backward(self, wavefront):
        """Propagate a wavefront backwards through the DM.
        
        Parameters
        ----------
        wavefront : Wavefront
            The incoming wavefront.
        
        Returns
        -------
        Wavefront
            The reflected wavefront.
        """
        wf = wavefront.copy()
        wf.electric_field *= np.exp(
                                 -2j*self.surfc*wavefront.wavenumber
                                   )
        return wf

    def phase_for(self, wavelength):
        """Get the phase that is added to a wavefront with a specified 
        wavelength.
        
        Parameters
        ----------
        wavelength : float
            The wavelength at which to calculate the phase deformation.
        
        Returns
        -------
        result : hci.Field
            The calculated phase deformation.
        """
        return 2*self.surfc*2*np.pi/wavelength


class TransmMask(hci.OpticalElement):
    """A `hci.OpticalElement` object representing a focal-plane mask.
    
    Parameters
    ----------
    mask_transm : ndarray
        Mask transmission function.  
        
    """
    def __init__(self, mask_transm):
        self.mask_transm = mask_transm

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
        the aperture of the deformable mirror (DM).
    acttr_filnm : str
        Fits filename for a single actuator influence function (2darray).  
        The same function is applied to all actuators. 
    save_filnm : str
        Filename for storing/loading array of influence functions in NPZ 
        format.  
    dm_diamtr: int
        Physical diameter of the DM pupil.  
    acttr_spacing : float
        Physical spacing between adjacent actuators on the DM.
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
    acttr : 2darray, @property
        Single-actuator influence function extracted from acttr_filnm
    acttr_grid : hci.grid, @property
        Automatically generated grid for a single actuator. Only used for 
        resampling/interpolation.
    modes : spr.csr_matrix, @property
        Peak-normalized DM influence functions.   
    
    Notes
    -----
    - Future plan:  Expand to use a separate influence function for each 
      actuator.  
 
    """
    def __init__(self, pupil_grid=hci.make_pupil_grid(481)
                 , acttr_filnm=os.path.join(os.path.expanduser('~'), 'sleec_in'
                                            , 'superinf-20180730-10x10.fits')
                 , save_filnm=os.path.join(os.path.expanduser('~')
                                           , 'sleec_out', 'inflnc_ftns.pkl')
                 , dm_diamtr=1.0, acttr_spacing=1.0/48
                 , acttr_pxl_spacing=10):
        self.pupil_grid = pupil_grid
        self.acttr_filnm = acttr_filnm
        self.save_filnm = save_filnm 
        self.dm_diamtr = dm_diamtr
        self.acttr_spacing = acttr_spacing
        self.acttr_pxl_spacing = acttr_pxl_spacing
    
    def auto_acttr_positns(self):
        """
        Generate an automatic, uniform grid of actuator positions.  

        Returns
        -------
        A `hci.Grid` of actuator positions.

        Notes
        -----
        - For use with resampling/interpolation.
 
        """  
        extent = self.dm_diamtr
        return hci.make_uniform_grid(self.dm_diamtr/self.acttr_spacing
                                     , [extent]*2)

    def make_intrpltr(self, field, grid=None, fill_value=np.nan
                      , method='nearest'):
        """Make an interpolator for a separated grid.  This method is 
        basically the hci.make_linear_interpolator_separated method, with the 
        added `method` option so one can choose linear or nearest interpolation.

        Parameters
        ----------
        field : hci.Field or ndarray
            The field to interpolate
        grid : hci.Grid or None
            The grid of the field.  If it is given, the grid of `field` is 
            replaced by this grid.  
        fill_value : scalar
            The value to use for points outside of the domain of the input 
            field.  If this is None, the values outside the domain are 
            extrapolated. 
        method : {'linear', 'nearest'}
            The method of interpolation to perform.  
 
        Returns
        -------
        hci.Field generator
            The interpolator, as a hci.Field generator.  The grid on which 
            this hci.Field generator will be evaluated, does not have to be 
            separated.  

        Notes
        -----
        One should only use 'nearest' option when `self.acttr_grid` 
        coordinates coincide with a subset of the `evaluated_grid` coordinates. 

        """
        if grid is None:
            grid = field.grid
        else:
            field = hci.Field(field, grid)
        axes_reversed = np.array(grid.separated_coords)
        interp = scint.RegularGridInterpolator(
                     points=axes_reversed, values=field.shaped, method=method
                     , bounds_error=False, fill_value=fill_value
                                              )
        def interpolator(evaluated_grid):
            evaluated_coords = np.flip(np.array(evaluated_grid.coords), 0)
            res = interp(evaluated_coords.T)
            return hci.Field(res.ravel(), evaluated_grid)
        return interpolator

    def genrt_modes(self, x_tilt=0.0, y_tilt=0.0, z_tilt=0.0, save=True
                 , intrpltn='nearest'):
        """Create DM influence functions.  If `save`=True, save the influence 
        functions as CSR sparse matrices.  If not, return a hci.ModeBasis 
        object representing the array of influence functions. 

        Parameters
        ----------
        x_tilt : float
            Tilt of the DM around the x-axis in radians.
        y_tilt : float
            Tilt of the DM around the y-axis in radians.
        z_tilt : float
            Tilt of the DM around the z-axis in radians.
        save : bool
            True --> save the array of influence functions as a CSR sparse 
            matrix in npz file format.  
            False --> return a hci.ModeBasis object representing the 
            array of influence functions.
        intrpltn : 'nearest' or 'linear'
            Interpolation method   

        Returns
        -------
        result : csr_matrix, if `save` = False
            Sparse matrix of peak-normalized influence functions.  If 
            `save`=True, a dictionary of the CSR sparse matrix of influence 
            functions (key = 'modes') and relevant paramters are pickled to a 
            file (filename = `self.save_filnm`).   
    
        Notes
        -----
        Comparison of saving with NPZ vs. pickle:  
            In our own benchmark tests of 48x48 actuators and 481x481 pixel 
            pupil grid, pickle required 1.1 s to load (on average), and it 
            produced a file size of 131 MB. NPZ required 1.9 s to load, with 
            a file size of 44.8 MB.   
        
        """
        evaluated_grid = (self.pupil_grid.scaled(1.0/np.cos([y_tilt, x_tilt]))
                         .rotated(-z_tilt))
        self.acttr_positns = self.auto_acttr_positns()
        modes = spr.lil_matrix((self.acttr_positns.size
                                , evaluated_grid.size))
        if intrpltn == 'nearest':
            acttr_out = self.make_intrpltr(
                            field=self.acttr.ravel(), grid=self.acttr_grid
                            , fill_value=0.0, method='nearest'
                                     )
        elif intrpltn == 'linear':
            acttr_out = self.make_intrpltr(
                            field=self.acttr.ravel(), grid=self.acttr_grid
                            , fill_value=0.0, method='linear'
                                     )
        for i in range(self.acttr_positns.size):
            modes[i, ] = acttr_out(evaluated_grid
                               .shifted(-self.acttr_positns[i]))
        if save:
#            spr.save_npz(self.save_filnm, modes.tocsr())
            save_dict = {'x_tilt': x_tilt, 'y_tilt': y_tilt, 'z_tilt': z_tilt
                         , 'intrpltn': intrpltn
                         , 'evaluated_grid': evaluated_grid
                         , 'self.acttr_positns': self.acttr_positns
                         , 'self.acttr_grid': self.acttr_grid
                         , 'self.acttr': self.acttr
                         , 'modes': modes.tocsr()} 
            f = open(self.save_filnm, 'wb')
            pickle.dump(save_dict, f)
            f.close()
        else:
            return modes

    def make_surfc(self, actuatn_map):
        """Return the DM surface in meters.  
 
        Parameters
        ----------
        actuatn_map : ndarray
            Coefficient (in meters) applied to each actuator (basically the 
            voltage map x m/V).  Rows correspond to actuators (e.g. 4096 rows 
            for 48x48 actuators.  Columns correspond to pixels (e.g. 231361 
            columns for 481x481 pixels).  

        Returns
        -------
        result : hci.Field
            The DM surface in meters.
      
        """
        surfc = self.modes.transpose().dot(actuatn_map)
        return hci.Field(surfc, self.pupil_grid)      

    def load_modes(self):
        f = open(self.save_filnm, 'rb')
        dict_results = pickle.load(f)
        return dict_results['modes']

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
    def acttr_positns(self):
        return self._acttr_positns
 
    @acttr_positns.setter
    def acttr_positns(self, user_input):
        self._acttr_positns = user_input

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

    @property
    def modes(self):
        return self._modes

    @modes.setter
    def modes(self, user_input):
        self._modes = user_input


class OptclSysPlanesOut(hci.OpticalElement):
    """Same as hci.OpticalSystem but returns the propagated wavefronts from 
    all optical elements.

    """
    def __init__(self, optical_elements):
        self.optical_elements = optical_elements
        self.wavefronts = []

    def forward(self, wavefront):
        """
        Forward propagation.  
    
        Parameters
        ----------
        wavefront : hci.Wavefront
            Input wavefront.
    
        Returns
        -------
        wavefronts : list
            Propagated wavefront from every optical element.  
        
        """
        wf = wavefront.copy()
        for optical_element in self.optical_elements:
            wf = optical_element.forward(wf)
            self.wavefronts.append(wf)
        return self.wavefronts
    
    def backward(self, wavefront):
        wf = wavefront.copy()
        for optical_element in reversed(self.optical_elements):
            wf = optical_element.backward(wf)
            self.wavefronts.append(wf)
        return self.wavefronts
    
    def get_transformation_matrix_forward(self, wavelength=1):
        matrix = 1
        for optical_element in self.optical_elements:
            matrix = np.dot(
                         optical_element.get_transformation_matrix_forward(
                                             wavelength), matrix
                           )
        return matrix
    
    def get_transformation_matrix_backward(self, wavelength=1):
        matrix = 1

        for optical_element in reversed(self.optical_elements):
            matrix = np.dot(
                         optical_element.get_transformation_matrix_backward(
                                             wavelength), matrix
                           )
        return matrix
    
    @property
    def optical_elements(self):
        return self._optical_elements

    @optical_elements.setter
    def optical_elements(self, optical_elements):
        self._optical_elements = list(optical_elements)
            

class Coronagraph(OptclSysPlanesOut):
    """SLEEC coronagraph model. 

    Notes
    -----
    - Need optimizable attributes to be modifiable from one iteration to the 
      next.    
    - Use @property to define those attributes.  
    - Have a method that imports <sleec_inp.py> to initialize the coronagraph 
    - Can have a function (outside of this class) to manipulate the attributes.

    """
    def __init__(self, wvl=None, pupil_field=None, fpm_field=None
                 , lyot_field=None, dark_field=None, dm1_vmap=None
                 , dm2_vmap=None, inflnc_ftn=None, dm_seprtn=None
                 , oversmpl=None, optical_elements=[]):
        super().__init__(optical_elements)
        self.wvl = wvl
        self.pupil_field = pupil_field
        self.fpm_field = fpm_field
        self.lyot_field = lyot_field
        self.dark_field = dark_field
        self.dm1_vmap = dm1_vmap
        self.dm2_vmap = dm2_vmap
        self.inflnc_ftn = inflnc_ftn
        self.dm_seprtn = dm_seprtn
        self.oversmpl = oversmpl

        self.optical_elements = [None]*10
        self.names = ['0: pupil', '1: dm1', '2: on_dm2', '3: dm2'
                      , '4: relay_pupil', '5: on_fpm', '6: fpm', '7: on_lyot'
                      , '8: lyot', '9: camera'] 

    def load_elements(self):
        self.optical_elements[0] = self.pupil
        self.optical_elements[1] = self.dm1
        self.optical_elements[2] = self.on_dm2
        self.optical_elements[3] = self.dm2
        self.optical_elements[4] = self.relay_pupil
        self.optical_elements[5] = self.on_fpm
        self.optical_elements[6] = self.fpm
        self.optical_elements[7] = self.on_lyot
        self.optical_elements[8] = self.lyot
        self.optical_elements[9] = self.camera

    # Pupil
    @property
    def pupil(self):
        pupil = TransmMask(self.pupil_field)
        return pupil
    # DM1
    @property
    def dm1(self):
        dm1 = DefrmblMirror(self.inflnc_ftn.make_surfc(self.dm1_vmap))
        return dm1
    # Propagator to DM2
    @property
    def on_dm2(self):
        prop_to_dm2 = hci.FresnelPropagator(self.pupil_field.grid
                                            , self.dm_seprtn)
        return prop_to_dm2
    # DM2
    @property
    def dm2(self):
        dm2 = DefrmblMirror(self.inflnc_ftn.make_surfc(self.dm2_vmap))
        return dm2
    # Propagator to Pupil
    @property
    def relay_pupil(self):
        prop_to_pupil = FresnelBack(self.pupil_field.grid, self.dm_seprtn)
        return prop_to_pupil
    # Propagator to Focal Plane
    @property
    def on_fpm(self):
        prop_to_focal = hci.FraunhoferPropagator(self.pupil_field.grid
                                                 , self.fpm_field.grid)
        return prop_to_focal
    # Focal-Plane Mask
    @property
    def fpm(self):
        fpm = TransmMask(self.fpm_field)
        return fpm
    # Propagator to Lyot Plane
    @property
    def on_lyot(self):
        prop_to_lyot = FraunhoferBack(self.pupil_field.grid
                                      , self.fpm_field.grid)
        return prop_to_lyot
    # Lyot-Plane Mask
    @property
    def lyot(self):
        lyot_mask = TransmMask(self.lyot_field)
        return lyot_mask
    # Propagator to Camera
    @property
    def camera(self):
        prop_to_camera = hci.FraunhoferPropagator(self.pupil_field.grid
                                                  , self.fpm_field.grid)
        return prop_to_camera
      

def demo_coronagraph():
    import matplotlib.pyplot as plt
    corngrph = Coronagraph()
    corngrph.load_args()
    corngrph.load_elements()
    output = corngrph(hci.Wavefront(corngrph.pupil_field))
    print(corngrph.names)
    hci.imshow_field(output[9].intensity, cmap='inferno')
    plt.colorbar()
    plt.show()

if __name__ == '__main__':
    demo_coronagraph()
