import numpy as np
import hcipy as hci


def fpm_dot(fpm_wid, corngrph, num_airy=28, field_stop=True):
    fpm_opacity = hci.circular_aperture(fpm_wid)
    if field_stop:
        fpm_transm = (
            (1.0 - fpm_opacity(corngrph.fpm_field.grid))
            *hci.circular_aperture(
                num_airy*corngrph.f_num*corngrph.wvl
                )(corngrph.fpm_field.grid)
        )
    else:
        fpm_transm = 1.0 - fpm_opacity(corngrph.fpm_field.grid)
    update_fpm(fpm_transm, corngrph)
    return corngrph


def lyot_circ_clear(lyot_wid, corngrph):
    lyot_transm = hci.circular_aperture(lyot_wid)(corngrph.lyot_field.grid)
    update_lyot(lyot_transm, corngrph)
    return corngrph


def update_pupil(pupil_field, corngrph):
    corngrph.pupil_field = pupil_field
    corngrph.optical_elements[0] = corngrph.pupil
    return corngrph


def update_dm1_surfc(surfc, corngrph):
    """Return the coronagraph with DM1's surface figure updated.

    Parameters
    ----------
    surfc : hcipy.Field
        DM surface figure
    corngrph : sleec_opt.Coronagraph
        The coronagraph to be updated

    Returns
    ------
    return : sleec_opt.Coronagraph
        The coronagraph with its DM1 updated

    """
    corngrph.dm1_surfc = surfc
    corngrph.optical_elements[1] = corngrph.dm1
    return corngrph


def update_dm2_surfc(surfc, corngrph):
    """Return the coronagraph with DM2's surface figure updated.

    Parameters
    ----------
    surfc : hcipy.Field
        DM surface figure
    corngrph : sleec_opt.Coronagraph
        The coronagraph to be updated

    Returns
    ------
    return : sleec_opt.Coronagraph
        The coronagraph with its DM1 updated

    """
    corngrph.dm2_surfc = surfc
    corngrph.optical_elements[3] = corngrph.dm2
    return corngrph


def update_dm1_vmap(vmap, corngrph):
    corngrph.dm1_vmap = vmap
    corngrph.optical_elements[1] = corngrph.dm1
    return corngrph


def update_dm2_vmap(vmap, corngrph):
    corngrph.dm2_vmap = vmap
    corngrph.optical_elements[3] = corngrph.dm2
    return corngrph


def update_fpm(fpm_field, corngrph):
    corngrph.fpm_field = fpm_field
    corngrph.optical_elements[6] = corngrph.fpm
    return corngrph


def update_lyot(lyot_field, corngrph):
    corngrph.lyot_field = lyot_field
    corngrph.optical_elements[8] = corngrph.lyot
    return corngrph


#def sphrcl_vmap(shape, focal_length, diamtr):
#    grid = hci.make_uniform_grid(dims=shape, extent=1)
#    pgrid = grid.as_('polar')
#    roc = focal_length*2.0
#    dist_from_ctr = np.sqrt(roc**2 - (diamtr/2.0)**2)
#    height = np.sqrt(dist_from_ctr**2 + pgrid.r**2) - dist_from_ctr
#    return height 


def _demo_dmsurfc():
    t0 = time.time()   
    coro = ini_corngrph(dm_surfc_args())
    t1 = time.time()
    print("initialization took {} s".format(t1-t0))
    output = coro(hci.Wavefront(coro.pupil_field, coro.wvl))
    print("propagaton took {} s".format(time.time()-t1))
    subplot_seq = 1
    plt.figure(figsize=(15, 5))
    for key in coro.names:
        print(key)
        wf = output[key]
        print(wf)
        plt.subplot(2, len(coro.names)//2, subplot_seq)
        hci.imshow_field(wf.intensity, cmap=turbo)
        plt.colorbar()
        subplot_seq += 1
    plt.show()
    plt.close()
    t0 = time.time()   
    
    
def _test_update_dm1(): 
    import matplotlib.pyplot as plt
    coro = ini_corngrph(dm_surfc_args())
    print(coro.names)
    wf0 = hci.Wavefront(electric_field=coro.pupil_field
                                    , wavelength=coro.wvl)
    output0 = coro(wf0)
    subplot_seq = 1
    plt.figure(figsize=(15, 5))
    for key in coro.names:
        print(key)
        wf = output0[key]
        print(wf)
        plt.subplot(2, len(coro.names)//2, subplot_seq)
        hci.imshow_field(wf.intensity, cmap=turbo)
        plt.colorbar()
        subplot_seq += 1
    plt.show()
    plt.close()
    surfc = hci.Field((coro.wvl/10)
                *np.cos(2*np.pi*(5/coro.pupil_diamtr)*coro.pupil_field.grid.x)
                , coro.pupil_field.grid) 
    update_dm1_surfc(surfc, coro)
    output1 = coro(wf0)
    subplot_seq = 1
    plt.figure(figsize=(15, 5))
    for key in coro.names:
        print(key)
        wf = output1[key]
        print(wf)
        plt.subplot(2, len(coro.names)//2, subplot_seq)
        hci.imshow_field(wf.intensity, cmap=turbo)
        plt.colorbar()
        subplot_seq += 1
    plt.show()
    plt.close()


def _test_update_dm2(): 
    import matplotlib.pyplot as plt
#    coro = ini_corngrph(dm_surfc_args())
    coro = ini_corngrph(default_args())
    print(coro.names)
    wf0 = hci.Wavefront(electric_field=coro.pupil_field
                                    , wavelength=coro.wvl)
    output0 = coro(wf0)
    surfc = hci.Field((coro.wvl/10)
                *np.cos(2*np.pi*(5/coro.pupil_diamtr)*coro.pupil_field.grid.y)
                , coro.pupil_field.grid) 
    update_dm2_surfc(surfc, coro)
    output1 = coro(wf0)
    subplot_seq = 1
    plt.figure(figsize=(15, 5))
    for key in coro.names:
        print(key)
        wf = output1[key]
        print(wf)
        plt.subplot(2, len(coro.names)//2, subplot_seq)
        hci.imshow_field(wf.intensity, cmap=turbo)
        plt.colorbar()
        subplot_seq += 1
    plt.show()
    plt.close()
