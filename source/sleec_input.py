"""Module to provide input arguments for coro_sleec.Coronagraph."""

import os
import numpy as np
import hcipy.hcipy as hci
import coro_sleec as cor

pars = {}


# wvl
pars['wvl'] = 1.0


# Generate `pupil_field`
num_acttrs_across = 48  # Number of actuators across the DM aperture.
acttr_pxl_spacing = 10  # Pixel sampling of the super-influence function, 
    # normally 10 pixels between adjacent actuators.  
pupil_dia_pxl = num_acttrs_across*acttr_pxl_spacing + 1  # Number of pixels 
    # across the pupil.
pupil_grid = hci.make_pupil_grid(pupil_dia_pxl)
pupil_transm = hci.circular_aperture(1)(pupil_grid)
pars['pupil_field'] = cor.TransmMask(pupil_grid, pupil_transm)


# Generate `inflnc_ftn'
acttr_filnm = os.path.join('..', '..', 'sleec_in'
                           , 'superinf-20180730-10x10.fits')  # Fits filename 
    # for a single-actuator influence function (2darray).  The same function 
    # is applied to all actuators.  
save_filnm = os.path.join('..', '..', 'sleec_out', 'inflnc_ftns.pkl')  # File 
    # name for a storing/loading the array of influence functions in NPZ 
    # format.  
dm_diamtr = 1.0  # Physical diameter of the DM
acttr_spacing = dm_diamtr/num_acttrs_across
inflnc_ftn = cor.InflncFtn(pupil_grid, acttr_filnm, save_filnm, dm_diamtr
                           , acttr_spacing, acttr_pxl_spacing)
try:
    inflnc_ftn.modes = inflnc_ftn.load_modes()
except FileNotFoundError:
    print("Influence-functions pickle file not found.")
    print("Generate influence functions on the fly.")
    inflnc_ftn.modes = inflnc_ftn.genrt_modes()
pars['inflnc_ftn'] = inflnc_ftn


#  `dm1_vmap`
pars['dm1_vmap'] = np.zeros(num_acttrs_across**2)


#  `dm2_vmap`
pars['dm2_vmap'] = np.zeros(num_acttrs_across**2)


#  `dm_seprtn`
pars['dm_seprtn'] = 1.0/wvl


# `oversmpl`
oversmpl = 8 # Pixel oversampling factor  
pars['oversmpl'] = oversmpl


# `fpm_field`
focal_grid = hci.make_focal_grid(pupil_grid, oversmpl, num_airy)
opacity = hci.circular_aperture(fpm_wid)
fpm_transm = 1.0 - fpm_opacity(focal_grid)  # Simple mask
pars['fpm_field'] = TransmMask(focal_grid, fpm_transm)

