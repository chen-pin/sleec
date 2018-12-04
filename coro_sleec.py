"""Module for modeling a coronagraph"""


i#from hcipy import hcipy as hci
import numpy as np
import matplotlib.pyplot as plt


class FocalPlaneMask(hci.OpticalElement):
    """A <hcipy> :class:'OpticalElement' object representing a focal-plane mask.
    
    Parmeters
    ---------
    input_grid : Grid
        The grid on which the incoming wavefront is defined.
    transm_field : ndarray
        Mask transmission function.  
    """
    def __init__(self, input_grid, transm_field):
        self.input_grid = input_grid
        self.transm_field = transm_field
    def forward(self, wavefront):
        print("wavefront\n{}\n".format(wavefront))
        wavefront.electric_field *= self.transm_field
        return wavefront
    def backward(self, wavefront):
        wavefront.electric_field *= self.transm_field
        return wavefront

class LyotMask(hci.OpticalElement):
    """
    """


class SleecCoronagraph(hci.OpticalElement):
    """A Lyot coronagraph.
    Parameters
    ----------
    input_grid : hci.Grid
        The grid on which the incoming wavefront is defined.
    focal_plane_mask : hci.OpticalElement
	:class:'OpticalElement' specifying the (complex) transmission of 
        the focal-plane mask. 
    lyot_stop : Field or OpticalElement or None
        :class:'OpticalElement' specifying the (complex) transmission of 
        the Lyot stop. 
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


def _test():  
    pupil_grid = hci.make_pupil_grid(128)
    focal_grid = hci.make_focal_grid(pupil_grid, 8, 32)
    prop = hci.FraunhoferPropagator(pupil_grid, focal_grid)
    aperture = hci.circular_aperture(6.0)
    fpm_transm = 1.0 - aperture(focal_grid)
    fpm = FocalPlaneMask(focal_grid, fpm_transm)
    pupil = hci.Wavefront(hci.circular_aperture(1)(pupil_grid))
    coro = SleecCoronagraph(pupil_grid, fpm)
    lyot_plane = coro(pupil)
    hci.imshow_field(lyot_plane.intensity)
    plt.colorbar()
    plt.show()
    lyot_stop = hci.Apodizer(hci.circular_aperture(254.0/256)(pupil_grid))
    prop = hci.FraunhoferPropagator(pupil_grid, focal_grid)
    img = prop(lyot_stop(lyot_plane))
    img_ref = prop(pupil)
    hci.imshow_field(lyot_stop(lyot_plane).intensity)
    plt.colorbar()
    plt.show()
    hci.imshow_field(np.log10(img.intensity/img_ref.intensity.max()), vmin=-12)
    plt.colorbar()
    plt.show()
#    dm = hci.make_xinetics_influence_functions(pupil_grid, 48, 1.0/47.)
#    plt.imshow(dm[2304//2+24].reshape(128,128))
#    plt.show()

