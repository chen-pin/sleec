"""<sleec_alter>:  Module for applying changes to (evolving the state of) 
a coronagraph."""

import os
import math
import pickle as pkl
import numpy as np
import hcipy as hci


def pupil_slepian_r(coefs, corngrph, file_exists, filename, m, n2d, k, L=170
                  , method='SE', orthgnlty=0, num_dense_pts=4096):
    """
    Return a coronagraph with pupil amplitude/apodization contructed by a 
    superposition of radial Slepian functions.  

    Parameters
    ----------
    coefs : array_like
        Basis function coefficients for the superposition
    corngrph : `sleec_optics.Coronagraph`
        The coronagraph object to be altered
    file_exists : boolean
        If `True`, use an existing file containing the mode basis of Slepians. 
        Otherwise, compute the mode basis and save it using `filename`.  
    filename : str
        If `file_exitsts` = `True`, it is the name of the existing mode-basis 
        file.  Othewise, it is the name of the save file.  
    m : int | None (if `file_exists` = `True`)
        The fixed order of the eigenvalue problem.  It is analogous to the 
        subscript n for Bessel functions J_n.
    n2d : float | None (if `file_exists` = `True`)
        The two-dimensional Shannon number.  For a Fraunhofer propagation, 
        n2d = (pi*N/4)^2, where N is the number of cycles across the 
        aperture. For example, if we are solving for a pupil-apodizing 
        function to limit energy in a 3wvl/D (radius) zone, n2d = 5.6.   
    k : int | None (if `file_exists` = `True`)
        The number of functions to output.  For example, k=3, the output 
        contains 3 Slepain functions g_1, g_2, g_3.  The g_k functions have 
        increasing frequencies with increasing k.  Examining a small number 
        examples suggest the number of zeros in the range [0, 1] equals k-1. 
        Note that k > n2d would not make much sense because the Shannon number  
        signifies the number of significant eigenfunctions (Slepians) that 
        can simultaneously be concentrated into the defined spatio-spectral 
        domain.
    L : int | None (if `file_exists` = `True`)
        The maximum degree in the Jacobi/Bessel expansion (default:  2*n2d),  
        or else the maximum polynomial degree for which GL wil be exact. 
        Simons states that increasing `L` is almost always a good idea,
        presumably because it increases accuracy.  However, P. Chen found that 
        computations can fail for high `L` values.
    method : str | None (if `file_exists` = `True`)
        'DV' compute by the method of De Villiers
        'SE' by Slepian "extension" [the default method]
        'GL' by direct Gauss-Legendre integration
    orthgnlty : int|None (if `file_exists` = `True`)
        1 Scales the solution for weightless orthogonality
        0 Scales the solution for orthogonality with weight r
    num_dense_pts : int | None (if `file_exists` = `True`)
        Number of densly sampled points in the Slepian computation
    
    Returns
    -------
    return : `hci.ModeBasis`
        A list of the user specified Slepian functions evaluated on the input 
        grid. 

    Notes
    -----
    An `hcipy.ModeBasis` object is a list of `hcipy.Field` objects.  

    """
    if file_exists: 
        with open(filename, 'rb') as f:
            slepians = pkl.load(f)['slepians']
    else:
        cwd = os.getcwd()
        os.chdir('slepian')
        from slepian.gen_slepians import gen_slepians_r as gsl
        slepians = gsl(m, n2d, k, L, method, orthgnlty, corngrph.pupil_diamtr
                       , corngrph.pupil_field.grid, num_dense_pts)
        os.chdir(cwd)
        with  open(filename, 'wb') as f:
            output = {}
            output['m'] = m
            output['n2d'] = n2d
            output['k'] = k
            output['L'] = L
            output['orthgnlty'] = orthgnlty
            output['D'] = corngrph.pupil_diamtr
            output['slepians'] = slepians
            pkl.dump(output, f)
    pupil_field = slepians.linear_combination(coefs)
    return update_pupil(pupil_field, corngrph) 
        

def fpm_dot(fpm_dia, corngrph, num_airy=28, field_stop=True):
    fpm_opacity = hci.circular_aperture(fpm_dia)
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


def lyot_circ_clear(lyot_dia, corngrph):
    lyot_transm = hci.circular_aperture(lyot_dia)(corngrph.lyot_field.grid)
    return update_lyot(lyot_transm, corngrph)


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


def _test_pupil_slepian_r():
    import yaml as yml
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    from turbo_colormap import turbo_colormap_data
    turbo = ListedColormap(turbo_colormap_data)
    import sleec_input as sli

    with open('sleec_pars.yml', 'r') as f:
        pars = yml.full_load(f)
    print('Loaded parameters from .yml file')
    corngrph = sli.import_pars()
    print('Coronagraph initialized')
    iwa = pars['fpm_dia']/2.0
    print('IWA = {}'.format(iwa))
    n2d = (np.pi*iwa/2.0)**2
    print('N2D = {}'.format(n2d))
    m = 0
    k = int(math.ceil(n2d))
#    k = 1
    print('k = {}'.format(k))
    coefs = np.random.random(k)
#    coefs = [1.0]
    print('Coefficients = {}'.format(coefs)) 
    file_exists = False
    filename = '../../sleec_out/slepians.pkl'
    corngrph = pupil_slepian_r(coefs, corngrph, file_exists, filename, m, n2d
                                 , k, L=170, method='SE', orthgnlty=0
                                 , num_dense_pts=4096)
    aperture = ( hci.circular_aperture(corngrph.pupil_diamtr)
                                    (corngrph.pupil_field.grid) )
    input_field= (
        hci.Field(np.ones_like(corngrph.pupil_field.grid.x)
                  , corngrph.pupil_field.grid)*aperture
                 )
    input_wavefront = hci.Wavefront(input_field, corngrph.wvl)
    output = corngrph(input_wavefront)
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


if __name__ == '__main__':
    _test_pupil_slepian_r()

