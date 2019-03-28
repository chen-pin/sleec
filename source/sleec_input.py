"""Module to provide input arguments for coro_sleec.Coronagraph."""

import os
import numpy as np
import hcipy.hcipy as hci
import coro_sleec as cor

def input_args():
    pars = {}
    
    # wvl
    wvl = 1.0
    pars['wvl'] = wvl
    
    # `pupil_field`
    num_acttrs_across = 48  # Number of actuators across the DM aperture.
    acttr_pxl_spacing = 10  # Pixel sampling of the super-influence function, 
        # normally 10 pixels between adjacent actuators.  
    pupil_dia_pxl = num_acttrs_across*acttr_pxl_spacing + 1  # Number of pixels 
        # across the pupil.
    pupil_grid = hci.make_pupil_grid(pupil_dia_pxl)
    pupil_field = hci.circular_aperture(1)(pupil_grid)
    pars['pupil_field'] = pupil_field
    
    # `inflnc_ftn'
    acttr_filnm = os.path.join('..', '..', 'sleec_in'
                               , 'superinf-20180730-10x10.fits')  # Fits 
                                                                  # filename 
        # for a single-actuator influence function (2darray).  The same function 
        # is applied to all actuators.  
    save_filnm = os.path.join('..', '..', 'sleec_out', 'inflnc_ftns.pkl')  
        # File 
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
    oversmpl = 8  # Pixel oversampling factor  
    pars['oversmpl'] = oversmpl
    
    # `fpm_field`
    num_airy = 24  
    fpm_wid = 6
    focal_grid = hci.make_focal_grid(pupil_grid, oversmpl, num_airy)
    fpm_opacity = hci.circular_aperture(fpm_wid)
    fpm_transm = 1.0 - fpm_opacity(focal_grid)  # Simple mask
    pars['fpm_field'] = fpm_transm
    
    # `lyot_field`
    lyot_dia = 0.9
    lyot_transm = hci.circular_aperture(lyot_dia)(pupil_grid)
    pars['lyot_field'] = lyot_transm
    
    # `dark_field`
    iwa = 3
    owa = 16
    pars['dark_field'] = ( hci.circular_aperture(owa)(focal_grid) 
                         - hci.circular_aperture(iwa)(focal_grid) )
    
    return pars


def ini_corngrph():
    # Initialize a coronagraph
    corngrph = cor.Coronagraph()
    args = input_args()
    corngrph.wvl = args['wvl']
    corngrph.pupil_field = args['pupil_field']
    corngrph.fpm_field = args['fpm_field']
    corngrph.lyot_field = args['lyot_field']
    corngrph.dark_field = args['dark_field']
    corngrph.dm1_vmap = args['dm1_vmap']
    corngrph.dm2_vmap = args['dm2_vmap']
    corngrph.inflnc_ftn = args['inflnc_ftn']
    corngrph.dm_seprtn = args['dm_seprtn']
    corngrph.oversmpl = args['oversmpl']
    corngrph.load_elements()
    return corngrph
