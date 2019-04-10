"""Module to provide input arguments for coro_sleec.Coronagraph."""

import os
import numpy as np
# import hcipy as hci
import hcipy_local.hcipy as hci
import sleec_optics as sop
import time

def input_args():
    pars = {}
    num_acttrs_across = 48  # Number of actuators across the DM aperture.
    acttr_pxl_spacing = 10  # Pixel sampling of the super-influence function, 
        # normally 10 pixels between adjacent actuators.  
    num_pxl_across = num_acttrs_across*acttr_pxl_spacing + 1  # Number of pixels 
        # across the pupil.
    acttr_spacing = 1E-3  # Actuator pitch
    
    # `wvl`
    wvl = 600E-9
    pars['wvl'] = wvl
    
    # `pupil_diamtr`
    pupil_diamtr = num_acttrs_across*acttr_spacing  # Physical pupil diameter
    pars['pupil_diamtr'] = pupil_diamtr
    pupil_grid = hci.make_pupil_grid(dims=num_pxl_across, diameter=pupil_diamtr)
    
    # `pupil_field` and
    pupil_field = hci.circular_aperture(diameter=pupil_diamtr)(pupil_grid)
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
    dm_diamtr = pupil_diamtr  # Physical diameter of the DM
    acttr_spacing = dm_diamtr/num_acttrs_across
    inflnc_ftn = sop.InflncFtn(pupil_grid, acttr_filnm, save_filnm, dm_diamtr
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
    pars['dm_seprtn'] = 1.0
    
    # `oversmpl`
    oversmpl = 8  # Pixel oversampling factor  
    pars['oversmpl'] = oversmpl
    
    # `f_num`    
    f_num = 30.0  # f/#
    pars['f_num'] = f_num

    # `fpm_field`
    num_airy = 28  # Outer extent in terms of Airy rings
    fpm_wid = 6*f_num*pars['wvl']  # FPM width/diameter
    focal_grid = hci.make_focal_grid(pupil_grid=pupil_grid, q=1
                                     , num_airy=num_airy
                                     , focal_length=pars['f_num']*pupil_diamtr
                                     , wavelength=pars['wvl'])
    fpm_opacity = hci.circular_aperture(fpm_wid)
    fpm_transm = 1.0 - fpm_opacity(focal_grid)  # Simple mask
    pars['fpm_field'] = fpm_transm
    
    # `lyot_field`
    lyot_dia = 0.9*pupil_diamtr
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
    corngrph = sop.Coronagraph()
    args = input_args()
    corngrph.wvl = args['wvl']
    corngrph.pupil_diamtr = args['pupil_diamtr']
    corngrph.pupil_field = args['pupil_field']
    corngrph.fpm_field = args['fpm_field']
    corngrph.lyot_field = args['lyot_field']
    corngrph.dark_field = args['dark_field']
    corngrph.dm1_vmap = args['dm1_vmap']
    corngrph.dm2_vmap = args['dm2_vmap']
    corngrph.inflnc_ftn = args['inflnc_ftn']
    corngrph.dm_seprtn = args['dm_seprtn']
    corngrph.oversmpl = args['oversmpl']
    corngrph.f_num = args['f_num']
    corngrph.load_elements()
    return corngrph


def pupil_lens(focal_length, diamtr, corngrph):
    pgrid = corngrph.pupil_field.grid.as_('polar')
    roc = focal_length*2.0
    dist_from_ctr = np.sqrt(roc**2 - (diamtr/2.0)**2)
    height = np.sqrt(dist_from_ctr**2 + pgrid.r**2) - dist_from_ctr
    corngrph.pupil_field = hci.Field(arr=height, grid=corngrph.pupil_field.grid)
    return corngrph 


def update_fpm(fpm_wid, corngrph):
    fpm_opacity = hci.circular_aperture(fpm_wid)
    fpm_transm = 1.0 - fpm_opacity(corngrph.fpm_field.grid)
    corngrph.fpm_field= fpm_transm
    corngrph.optical_elements[6] = corngrph.fpm
    return corngrph


def thrghpt_on_axis(prop_output):
    pupil_sum = np.sum(prop_output[0].intensity)
    dm2_sum = np.sum(prop_output[3].intensity)
    return dm2_sum/pupil_sum


def thrghpt_off_axis(corngrph, r=3.0, theta=np.pi/4.0, ampltd=0.1/(2*np.pi)):
    cgrid = corngrph.pupil_field.grid
    cartsn_to_polar = hci.field.polar_grid._cartesian_to_polar
    pgrid = cartsn_to_polar(cgrid)
    value = ampltd*np.exp(-1j*2*np.pi*pgrid.r*r*np.cos(theta-pgrid.theta))
    field = corngrph.pupil_field*(1.0+hci.Field(value, cgrid))
    wavefront = hci.Wavefront(field)
    return corngrph(wavefront)


def demo():
    import matplotlib.pyplot as plt
    t0 = time.time()   
    coro = ini_corngrph()
    output = coro(hci.Wavefront(coro.pupil_field))
    print("ini_corngrph & propagation took {} s".format(time.time()-t0))
    plt.subplot(221)
    hci.imshow_field(coro.fpm_field, cmap='inferno')
    plt.colorbar()
    plt.subplot(222)
    hci.imshow_field(output[9].intensity, cmap='inferno')
    plt.colorbar()
    t0 = time.time()   
    update_fpm(1.0, coro)
    output = coro(hci.Wavefront(coro.pupil_field))
    print("update_fpm and propagation took {} s".format(time.time()-t0))
    plt.subplot(223)
    hci.imshow_field(coro.fpm_field, cmap='inferno')
    plt.colorbar()
    plt.subplot(224)
    hci.imshow_field(output[9].intensity, cmap='inferno')
    plt.colorbar()
    plt.show()


if __name__ == '__main__':
    demo()
