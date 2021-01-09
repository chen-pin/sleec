"""<sleec_input>.  Module to load parameters and generate input 
arguments to initialize a coronagraph."""

import os
import time
import numpy as np
import yaml as yml
import astropy.io.fits as fits
import hcipy as hci
import sleec_optics as sop
import sleec_alter as sla


def import_pars(filename='sleec_pars.yml'):
    """Parse <sleec_pars.yml> and return an initialized coronagraph.  
    
    Returns
    -------
    output : `sleec_optics.coronagraph`
        Coronagraph object prescribed by <sleec_pars.yml>

    Notes
    -----
    Currently, it sets up only a classic Lyot coronagraph, with a clear 
    circular pupil, dot FPM, and a clear circular Lyot stop.

    """
    with open(filename, 'r') as f:    
        pars = yml.full_load(f)
    print('The following parameters have been loaded from ' + f.name)  
    for key, value in pars.items():  
        print("{}:  {}".format(key, value))
    
 
    # PUPIL	
    pupil_grid = hci.make_pupil_grid(dims=pars['num_pxl_across'], 
                                     diameter=pars['pupil_dia']
                                              *pars['grid_oversize'])
    pupil_field = ( hci.circular_aperture(diameter=pars['pupil_dia'])
                                           (pupil_grid) )

    # DM1
    if pars['dm1_surfc'] == 'zero':
        dm1_surfc = hci.Field(np.zeros(pars['num_pxl_across']**2)
                                      , pupil_grid)
        dm1_vmap = None
        dm2_vmap = None
        inflnc_ftn = None
    elif pars['dm1_surfc'] == None:
        acttr_path = os.path.join('..', '..', 'sleec_in'
                                   , pars['acttr_filnm'])  
        acttr_save_path = os.path.join('..', '..', 'sleec_out'
                                  , pars['actter_save_filnm'])  
        dm_dia = pars['pupil_dia']  
        acttr_spacing = dm_dia/pars['num_acttrs_across']
        inflnc_ftn = sop.InflncFtn(pupil_grid, acttr_path
                                   , acttr_save_path
                                   , dm_dia, acttr_spacing
                                   , acttr_pxl_spacing)
        try:
            inflnc_ftn.modes = inflnc_ftn.load_modes()
        except FileNotFoundError:
            print("Influence-functions pickle file not found.")
            print("Generate influence functions on the fly.")
            inflnc_ftn.modes = inflnc_ftn.genrt_modes()
        if dm1_vmap != None:  
            fn = os.path.join('..', '..', 'sleec_in', pars['dm1_vmap'])
            dm1_vmap = fits.getdata(fn, ext=0)
        else:  
            dm1_vmap = np.zeros(num_acttrs_across**2)
    else:
        fn = os.path.join('..', '..', 'sleec_in', pars['dm1_surfc'])
        dm1_surfc = fits.getdata(fn, ext=0)
        dm1_vmap = None
        dm2_vmap = None
        inflnc_ftn = None

    # DM2
    if pars['dm2_surfc'] == 'zero':
        dm2_surfc = hci.Field(np.zeros(pars['num_pxl_across']**2)
                                      , pupil_grid)
        dm1_vmap = None
        dm2_vmap = None
        inflnc_ftn = None
    elif pars['dm2_surfc'] == None:
        acttr_path = os.path.join('..', '..', 'sleec_in'
                                   , pars['acttr_filnm'])  
        acttr_save_path = os.path.join('..', '..', 'sleec_out'
                                  , pars['actter_save_filnm'])  
        dm_dia = pars['pupil_dia']  
        acttr_spacing = dm_dia/pars['num_acttrs_across']
        inflnc_ftn = sop.InflncFtn(pupil_grid, acttr_path
                                   , acttr_save_path
                                   , dm_dia, acttr_spacing
                                   , acttr_pxl_spacing)
        try:
            inflnc_ftn.modes = inflnc_ftn.load_modes()
        except FileNotFoundError:
            print("Influence-functions pickle file not found.")
            print("Generate influence functions on the fly.")
            inflnc_ftn.modes = inflnc_ftn.genrt_modes()
        if dm2_vmap != None:  
            fn = os.path.join('..', '..', 'sleec_in', pars['dm2_vmap'])
            dm2_vmap = fits.getdata(fn, ext=0)
        else:  
            dm2_vmap = np.zeros(num_acttrs_across**2)
    else:
        fn = os.path.join('..', '..', 'sleec_in', pars['dm2_surfc'])
        dm2_surfc = fits.getdata(fn, ext=0)
        dm1_vmap = None
        dm2_vmap = None
        inflnc_ftn = None

    
    # FOCAL PLANE
    focal_grid = hci.make_focal_grid(
                                     q=pars['oversmpl']
                                     , num_airy=pars['num_airy']
                                     , spatial_resolution
                                       =pars['f_num']*pars['wvl']
                                     )
#    if pars['fpm_type'] == 'dot':
#        fpm_opacity = hci.circular_aperture(pars['fpm_dia']
#                                            *pars['f_num']*pars['wvl'])
#        fpm_field = ( (1.0 - fpm_opacity(focal_grid))
#                       *hci.circular_aperture(
#                                2*pars['num_airy']*pars['f_num']*pars['wvl']
#                                             )(focal_grid)
#                     )
    
    fpm_opacity = hci.circular_aperture(pars['fpm_dia']
                                        *pars['f_num']*pars['wvl'])
    fpm_field = ( (1.0 - fpm_opacity(focal_grid))
                   *hci.circular_aperture(
                            2*pars['num_airy']*pars['f_num']*pars['wvl']
                                         )(focal_grid)
                 )

    # LYOT PLANE
    lyot_field = ( hci.circular_aperture(pars['lyot_dia']*pars['pupil_dia'])
                                      (pupil_grid) )
    
    # INITIALIZE CORONAGRAPH  
    args = {}
    args['wvl'] = pars['wvl']
    args['pupil_dia'] = pars['pupil_dia']
    args['pupil_field'] = pupil_field  
    args['fpm_field'] = fpm_field  
    args['lyot_field'] = lyot_field
    args['dm1_surfc'] = dm1_surfc  
    args['dm2_surfc'] = dm2_surfc
    args['dm1_vmap'] = dm1_vmap  
    args['dm2_vmap'] = dm2_vmap  
    args['inflnc_ftn'] = inflnc_ftn
    args['dm_seprtn'] = pars['dm_seprtn']  
    args['f_num'] = pars['f_num'] 
    coronagraph = ini_corngrph(args)
 
    return coronagraph, pars


def ini_corngrph(args):
    """ Initialize a coronagraph.  

    Parameters
    ----------
    args : dictionary  
        Input arguments for the `sleec_optics.Coronagraph` class.  

    Returns
    -------
    return : `sleec_optics.Coronagraph` object

    """ 
    corngrph = sop.Coronagraph()
    corngrph.wvl = args['wvl']
    corngrph.pupil_diameter = args['pupil_dia']
    corngrph.pupil_field = args['pupil_field']
    corngrph.fpm_field = args['fpm_field']
    corngrph.lyot_field = args['lyot_field']
    corngrph.dm1_surfc= args['dm1_surfc']
    corngrph.dm2_surfc = args['dm2_surfc']
    corngrph.dm1_vmap = args['dm1_vmap']
    corngrph.dm2_vmap = args['dm2_vmap']
    corngrph.inflnc_ftn = args['inflnc_ftn']
    corngrph.dm_seprtn = args['dm_seprtn']
    corngrph.f_num = args['f_num']
    corngrph.load_elements()
    return corngrph


def _demo():
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from turbo_colormap import turbo_colormap_data
    turbo = ListedColormap(turbo_colormap_data)
    t0 = time.time()   
    corngrph, pars = import_pars()
    output = corngrph(hci.Wavefront(corngrph.pupil_field, corngrph.wvl))
    print("ini_corngrph & propagation took {} s".format(time.time()-t0))
    subplot_seq = 1
    plt.figure(figsize=(15, 5))
    for key in corngrph.names:
        print(key)
        wf = output[key]
        print(wf)
        plt.subplot(2, len(corngrph.names)//2, subplot_seq)
        hci.imshow_field(wf.intensity, cmap=turbo)
        plt.colorbar()
        subplot_seq += 1
    plt.show()
    plt.close()
    contrast_map = (output['camera'].intensity
                    /output['on_fpm'].intensity.max())
    grid_scale = pars['f_num']*pars['wvl']
    hci.imshow_field(np.log10(contrast_map), cmap=turbo
                     ,vmax=0.0, vmin=-12.0, grid_units=grid_scale)
    plt.title('LOG CONTRAST')
    plt.colorbar()
    plt.show()
    plt.close()


if __name__ == '__main__':
    _demo()
