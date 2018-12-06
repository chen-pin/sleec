"""Module for modeling a coronagraph"""

import os
import matplotlib.pyplot as plt
import pickle
import time
import numpy as np
from hcipy import hcipy as hci


class TransmMask(hci.OpticalElement):
    """A `hci.OpticalElement` object representing a focal-plane mask.
    
    Parameters
    ----------
    input_grid : hci.Grid
        The grid on which the incoming wavefront is defined.
    mask_transm : ndarray
        Mask transmission function.  
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


class SleecCoronagraph(hci.OpticalElement):
    """A Lyot coronagraph.
    
    Parameters
    ----------
    input_grid : hci.Grid
        The grid on which the incoming wavefront is defined.
    focal_plane_mask : hci.OpticalElement
	Complex transmission of the focal-plane mask. 
    lyot_stop : Field or OpticalElement or None
        Complex transmission of the Lyot stop. 
    """
    
    def __init__(self, input_grid, focal_plane_mask, lyot_stop=None):
        self.focal_plane_mask = focal_plane_mask
        self.lyot_stop = lyot_stop
        self.prop = hci.FraunhoferPropagator(input_grid, focal_plane_mask
                                                         .input_grid)
        
    def forward(self, wavefront):
        wf_foc = self.prop.forward(wavefront)  
        wf_foc.electric_field = (self.focal_plane_mask.forward(wf_foc)
                                     .electric_field)
        lyot = self.prop.backward(wf_foc)        
        if self.lyot_stop is not None:
            lyot = self.lyot_stop.forward(lyot)            
        return lyot
    
    def backward(self, wavefront):
        if self.lyot_stop is not None:
            wf = self.lyot_stop.backward(wavefront)
        else:
            wf = wavefront            
        wf_foc = self.prop.forward(wf)
        wf_foc.electric_field = (self.focal_plane_mask.backward(wf_foc)
                                     .electric_field)       
        pup = self.prop.backward(wf_foc)
        pup.electric_field[:] = wavefront.electric_field - wf_foc.electric_field        
        return pup


def _test_class():  
    pupil_grid = hci.make_pupil_grid(256)
    focal_grid = hci.make_focal_grid(pupil_grid, 8, 128)
    prop = hci.FraunhoferPropagator(pupil_grid, focal_grid)
#    # Gaussian mask function
#    x, y = focal_grid.coords
#    r = np.sqrt(x**2 + y**2)
#    fpm_transm = 1.0 - np.exp(-(np.log(2.0)/24.0**2)*r**2)
    # Binary mask function
    aperture = hci.circular_aperture(6.5)
    fpm_transm = 1.0 - aperture(focal_grid)
    fpm = TransmMask(focal_grid, fpm_transm)
    pupil = hci.Wavefront(hci.circular_aperture(1)(pupil_grid))
    img_ref = prop(pupil)
    hci.imshow_field(fpm(img_ref).amplitude)
    plt.colorbar()
    plt.show()
    coro = SleecCoronagraph(pupil_grid, fpm)
    lyot_plane = coro(pupil)
#    print(lyot_plane.electric_field)
#    plt.imshow(np.square(np.abs(lyot_plane.electric_field.reshape(128,128))), cmap='inferno')
    hci.imshow_field(lyot_plane.intensity)
    plt.colorbar()
    plt.show()
    lyot_stop = hci.Apodizer(hci.circular_aperture(254.0/256)(pupil_grid))
    prop = hci.FraunhoferPropagator(pupil_grid, focal_grid)
    img = prop(lyot_stop(lyot_plane))
    hci.imshow_field(lyot_stop(lyot_plane).intensity)
    plt.colorbar()
    plt.show()
    hci.imshow_field(np.log10(img.intensity/img_ref.intensity.max()), vmin=-12)
    plt.colorbar()
    plt.show()
#    dm = hci.make_xinetics_influence_functions(pupil_grid, 48, 1.0/47.)
#    hci.imshow_field(dm[2304//2+24])
#    plt.colorbar()
#    plt.show()


def _test_instance(dm1_array=np.zeros(48*48), dm2_array=np.zeros(48*48)
                   , pupil_dia_pxl=256, fpm_wid=6.0, lyot_dia=0.9
                   , oversmpl=8, num_airy=24, dm_dia_num_acttr=48
                   , dm_seprtn=1/600e-6):
    t0 = time.time()
    optical_elements = []
    # Make pupil
    pupil_grid = hci.make_pupil_grid(pupil_dia_pxl)
    focal_grid = hci.make_focal_grid(pupil_grid, oversmpl, num_airy)
    pupil_transm = hci.circular_aperture(1)(pupil_grid)
    pupil = TransmMask(pupil_grid, pupil_transm)
    optical_elements.append(pupil)
    # Influence Functions
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
#    plt.show()
    plt.savefig('dm2_int.jpg')
    plt.close()
    hci.imshow_field(camera_wf.phase, vmin=-1, vmax=1, cmap='gray')
    plt.colorbar()
#    plt.show()
    plt.savefig('dm2_ph.jpg')
    plt.close()
#    reshape_wid = int(np.sqrt(camera_wf.phase.shape[0]))
#    diff = (camera_wf.phase.reshape((reshape_wid, -1))[reshape_wid//2, :]
#           - camera_wf.phase.reshape((reshape_wid, -1))[:, reshape_wid//2])
#    plt.plot(diff)
#    plt.show()
