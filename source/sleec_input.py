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
    acttr_spacing = 1E-3  # Actuator pitch 
    num_pxl_across = num_acttrs_across*acttr_pxl_spacing + 1  # Number of pixels 
        # across the pupil.
    
    # `wvl`
    wvl = 600E-9
    pars['wvl'] = wvl
    
    # `pupil_diamtr`
    pupil_diamtr = num_acttrs_across*acttr_spacing  # Physical pupil diameter
    pars['pupil_diamtr'] = pupil_diamtr
    pupil_grid = hci.make_pupil_grid(dims=num_pxl_across, diameter=pupil_diamtr)
    
    # `pupil_field` 
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
    
    # `f_num`    
    f_num = 30.0  # f/#
    pars['f_num'] = f_num

    # `fpm_field`
    num_airy = 28  # Outer extent in terms of Airy rings
    fpm_wid = 6.0*f_num*pars['wvl']  # FPM width/diameter
    oversmpl = 8  # Pixel oversampling factor 
    focal_grid = hci.make_focal_grid(pupil_grid=pupil_grid, q=oversmpl
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
    iwa = 3.0*f_num*pars['wvl']
    owa = 16*f_num*pars['wvl']
    pars['dark_field'] = ( hci.circular_aperture(2.0*owa)(focal_grid) 
                         - hci.circular_aperture(2.0*iwa)(focal_grid) )
    
    return pars


def default_args():
    pars = {}
    num_acttrs_across = 48  # Number of actuators across the DM aperture.
    acttr_pxl_spacing = 10  # Pixel sampling of the super-influence function, 
        # normally 10 pixels between adjacent actuators. 
    acttr_spacing = 1E-3  # Actuator pitch 
    num_pxl_across = num_acttrs_across*acttr_pxl_spacing + 1  # Number of pixels 
        # across the pupil.
    
    # `wvl`
    wvl = 600E-9
    pars['wvl'] = wvl
    
    # `pupil_diamtr`
    pupil_diamtr = num_acttrs_across*acttr_spacing  # Physical pupil diameter
    pars['pupil_diamtr'] = pupil_diamtr
    pupil_grid = hci.make_pupil_grid(dims=num_pxl_across, diameter=pupil_diamtr)
    
    # `pupil_field` 
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
    
    # `f_num`    
    f_num = 30.0  # f/#
    pars['f_num'] = f_num

    # `fpm_field`
    num_airy = 28  # Outer extent in terms of Airy rings
    fpm_wid = 6.0*f_num*pars['wvl']  # FPM width/diameter
    oversmpl = 8  # Pixel oversampling factor 
    focal_grid = hci.make_focal_grid(pupil_grid=pupil_grid, q=oversmpl
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
    iwa = 3.0*f_num*pars['wvl']
    owa = 16*f_num*pars['wvl']
    pars['dark_field'] = ( hci.circular_aperture(2.0*owa)(focal_grid) 
                         - hci.circular_aperture(2.0*iwa)(focal_grid) )
    
    return pars


def ini_corngrph(args=default_args()):
    # Initialize a coronagraph
    corngrph = sop.Coronagraph()
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
    corngrph.f_num = args['f_num']
    corngrph.load_elements()
    return corngrph


def sphrcl_vmap(shape, focal_length, diamtr):
    grid = hci.make_uniform_grid(dims=shape, extent=1)
    pgrid = grid.as_('polar')
    roc = focal_length*2.0
    dist_from_ctr = np.sqrt(roc**2 - (diamtr/2.0)**2)
    height = np.sqrt(dist_from_ctr**2 + pgrid.r**2) - dist_from_ctr
    return height 


def update_fpm(fpm_wid, corngrph):
    fpm_opacity = hci.circular_aperture(fpm_wid)
    fpm_transm = 1.0 - fpm_opacity(corngrph.fpm_field.grid)
    corngrph.fpm_field= fpm_transm
    corngrph.optical_elements[6] = corngrph.fpm
    return corngrph


def fpm_thrghpt(prop_output, corngrph):
    fpm_in = np.sum(prop_output[5].intensity*corngrph.dark_field)
    fpm_out = np.sum(prop_output[6].intensity*corngrph.dark_field)
    return fpm_out/fpm_in


def lyot_thrghpt(corngrph):
    indx = np.where(np.absolute(corngrph.pupil.mask_transm) != 0.0)
    total = np.sum(np.absolute(corngrph.pupil.mask_transm)[indx]**2)
    sum_transm_ratio = np.sum(
            np.square(np.absolute(corngrph.lyot.mask_transm)[indx]
                      /np.absolute(corngrph.pupil.mask_transm)[indx]
                     )
                             )
    return sum_transm_ratio/total


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


def _test_fpm_thrghpt():
    """Test <fpm_thrghpt> function.  To get valid results, 
    fpm_coef must be > 2*iwa and owa <= 28 (outer extent of focal grid)"""

    import matplotlib.pyplot as plt
    from scipy import special as spl
    corngrph = ini_corngrph()
    # Update FPM
    f_num = 30.0
    fpm_coef = 7.0
    fpm_wid = fpm_coef*f_num*corngrph.wvl  # FPM width/diameter
    fpm_opacity = hci.circular_aperture(fpm_wid)
    fpm_transm = 1.0 - fpm_opacity(corngrph.fpm_field.grid)  # Simple mask
    corngrph.fpm_field = fpm_transm
    # Update Dark Zone  
    iwa_coef = 3.0  
    iwa = iwa_coef*f_num*corngrph.wvl
    owa_coef = 28.0
    owa = owa_coef*f_num*corngrph.wvl
    corngrph.dark_field = ( hci.circular_aperture(2.0*owa)(corngrph.fpm_field.grid) 
                         - hci.circular_aperture(2.0*iwa)(corngrph.fpm_field.grid) )
    print(corngrph.names)
    wf0 = hci.Wavefront(electric_field=corngrph.pupil_field
                                    , wavelength=corngrph.wvl)
    output0 = corngrph(wf0)
    predicted_thrghpt = ( (spl.jv(0, np.pi*fpm_coef/2.0)**2 
                         + spl.jv(1, np.pi*fpm_coef/2.0)**2 
                           - spl.jv(0, np.pi*owa_coef)**2 
                           - spl.jv(1, np.pi*owa_coef)**2)
                         /(spl.jv(0, np.pi*iwa_coef)**2 
                           + spl.jv(1, np.pi*iwa_coef)**2 
                           - spl.jv(0, np.pi*owa_coef)**2 
                           - spl.jv(1, np.pi*owa_coef)**2) )
    print("Predicted throughput = {}".format(predicted_thrghpt))
    print("Simulated throughput = {}".format(fpm_thrghpt(output0, corngrph)))
    plt.figure(figsize=(15, 5))
    plt.subplot(221)
    hci.imshow_field(output0[5].intensity, cmap='inferno')
    plt.colorbar()
    plt.subplot(222)
    hci.imshow_field(output0[6].intensity, cmap='inferno')
    plt.colorbar()
    plt.subplot(223)
    hci.imshow_field(corngrph.dark_field, cmap='inferno')
    plt.colorbar()
    plt.subplot(224)
    hci.imshow_field(corngrph.fpm.mask_transm, cmap='inferno')
    plt.colorbar()
    plt.show()


def _test_lyot_thrghpt():
    """Test <lyot_thrghpt> function.  To get valid results, 
    lyot_coef must be < 1"""

    import matplotlib.pyplot as plt
    corngrph = ini_corngrph()
    # Update FPM
    f_num = 30.0
    # Update Lyot Stop
    lyot_coef = 0.5
    lyot_dia = lyot_coef*corngrph.pupil_diamtr
    lyot_transm = hci.circular_aperture(lyot_dia)(corngrph.pupil_field.grid)
    corngrph.lyot_field = lyot_transm
    print(corngrph.names)
    wf0 = hci.Wavefront(electric_field=corngrph.pupil_field
                                    , wavelength=corngrph.wvl)
    output0 = corngrph(wf0)
    predicted_thrghpt = lyot_coef**2
    print("Predicted throughput = {}".format(predicted_thrghpt))
    print("Simulated throughput = {}".format(lyot_thrghpt(corngrph)))
    plt.figure(figsize=(15, 5))
    plt.subplot(121)
    hci.imshow_field(np.absolute(corngrph.pupil.mask_transm)**2, cmap='inferno')
    plt.colorbar()
    plt.subplot(122)
    hci.imshow_field(np.absolute(corngrph.lyot.mask_transm)**2, cmap='inferno')
    plt.colorbar()
    plt.show()


if __name__ == '__main__':
    _test_lyot_thrghpt()
