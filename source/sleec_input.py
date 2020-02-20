"""Module to provide input arguments for coro_sleec.Coronagraph."""

import os
import numpy as np
# import hcipy as hci
import matplotlib.pyplot as plt
#import hcipy_local.hcipy as hci
import hcipy as hci
import sleec_optics as sop
import time

def input_args(oversmpl=4, num_airy=28, field_stop=True):
    pars = {}
    num_acttrs_across = 48  # Number of actuators across the DM aperture.
    acttr_pxl_spacing = 10  # Pixel sampling of the super-influence function, 
        # normally 10 pixels between adjacent actuators. 
    acttr_spacing = 1E-3  # Actuator pitch 
    num_pxl_across = 2*(num_acttrs_across*acttr_pxl_spacing + 1)  # Number of pixels 
        # across the pupil.
    # `wvl`
    wvl = 600E-9
    pars['wvl'] = wvl
    
    # `pupil_diamtr`
    pupil_diamtr = num_acttrs_across*acttr_spacing  # Physical pupil diameter
    pars['pupil_diamtr'] = pupil_diamtr
    pupil_grid = hci.make_pupil_grid(dims=num_pxl_across, diameter=pupil_diamtr)
    
    # `pupil_field` 
    pupil_field = hci.circular_aperture(diameter=0.98*pupil_diamtr)(pupil_grid)
    pars['pupil_field'] = pupil_field

    # `inflnc_ftn'
    acttr_filnm = os.path.join('..', '..', 'sleec_in'
                               , 'superinf-20180730-10x10.fits')  # Fits 
                                                                  # filename 
        # for a single-actuator influence function (2darray).  The same function 
        # is applied to all actuators.  
    save_filnm = os.path.join('..', '..', 'sleec_out', 'inflnc_ftns_2x.pkl')  
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
    # print("Generate influence functions on the fly.")
    # inflnc_ftn.modes = inflnc_ftn.genrt_modes(save=True, intrpltn='linear')
    # print("Done")
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
#    num_airy = 28  # Outer extent in terms of Airy rings
    # fpm_wid = 6.0*f_num*pars['wvl']  # FPM width/diameter
    fpm_wid = 5.94098759e-05
#    oversmpl = 4  # Pixel oversampling factor 
#    focal_grid = hci.make_focal_grid(pupil_grid=pupil_grid, q=oversmpl
#                                     , num_airy=num_airy
#                                     , focal_length=pars['f_num']*pupil_diamtr
#                                     , wavelength=pars['wvl'])
    focal_grid = hci.make_focal_grid_from_pupil_grid(
                                     pupil_grid=pupil_grid
                                     , q=oversmpl
                                     , num_airy=num_airy
                                     , focal_length=pars['f_num']*pupil_diamtr
                                     , wavelength=pars['wvl']
                                                    )
    fpm_opacity = hci.circular_aperture(fpm_wid)
    if field_stop:
        fpm_transm = (1.0 - fpm_opacity(focal_grid))*hci.circular_aperture(num_airy*pars['f_num']*pars['wvl'])(focal_grid)
    else:
        fpm_transm = 1.0 - fpm_opacity(focal_grid)  # Simple mask
    pars['fpm_field'] = fpm_transm
    
    # `lyot_field`
    # lyot_dia = 0.9*pupil_diamtr
    lyot_dia = 2.80100555e-02
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
    pupil_field = hci.circular_aperture(diameter=0.98*pupil_diamtr)(pupil_grid)
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
    # fpm_wid = 6.0*f_num*pars['wvl']  # FPM width/diameter
    fpm_wid = 5.94098759e-05
    oversmpl = 4  # Pixel oversampling factor 
#    focal_grid = hci.make_focal_grid(pupil_grid=pupil_grid, q=oversmpl
#                                     , num_airy=num_airy
#                                     , focal_length=pars['f_num']*pupil_diamtr
#                                     , wavelength=pars['wvl'])
    focal_grid = hci.make_focal_grid_from_pupil_grid(
                                     pupil_grid=pupil_grid
                                     , q=oversmpl
                                     , num_airy=num_airy
                                     , focal_length=pars['f_num']*pupil_diamtr
                                     , wavelength=pars['wvl']
                                                    )
    fpm_opacity = hci.circular_aperture(fpm_wid)
    # fpm_transm = 1.0 - fpm_opacity(focal_grid)  # Simple mask
    fpm_transm = (1.0 - fpm_opacity(focal_grid))*hci.circular_aperture(num_airy*pars['f_num']*pars['wvl'])(focal_grid)
    pars['fpm_field'] = fpm_transm
    
    # `lyot_field`
    # lyot_dia = 0.9*pupil_diamtr
    lyot_dia = 2.80100555e-02
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


def update_fpm(fpm_wid, corngrph, num_airy=28, field_stop=True):
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
    corngrph.fpm_field= fpm_transm
    corngrph.optical_elements[6] = corngrph.fpm
    return corngrph


def update_lyot(lyot_wid, corngrph):
    lyot_transm = hci.circular_aperture(lyot_wid)(corngrph.pupil_field.grid)
    corngrph.lyot_field = lyot_transm
    corngrph.optical_elements[8] = corngrph.lyot
    return corngrph


def fpm_thrghpt(prop_output, corngrph):
    fpm_in = np.sum(prop_output['on_fpm'].intensity*corngrph.dark_field)
    fpm_out = np.sum(prop_output['fpm'].intensity*corngrph.dark_field)
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


# def thrghpt_off_axis(corngrph, ref_img, r, theta=np.radians(60)
#                      , ampltd=0.1/(2*np.pi)):
#     cgrid = corngrph.pupil_field.grid
#     cartsn_to_polar = hci.field.polar_grid._cartesian_to_polar
#     pgrid = cartsn_to_polar(cgrid)
#     value = ampltd*np.exp(-1j*2*np.pi*pgrid.r*r*np.cos(theta-pgrid.theta))
#     field = corngrph.pupil_field*(1.0+hci.Field(value, cgrid))
#     wavefront = hci.Wavefront(field)
#     img = corngrph(wavefront)['on_fpm'].intensity
#     return hci.get_strehl_from_focal(img, ref_img)


def strehl(prop_output, ref_intnsty): 
    return float(prop_output['on_fpm'].intensity.max()/ref_intnsty)


def mean_contrast(corngrph, prop_output, ref_intnsty):
    contrast_img = (prop_output['camera'].intensity
                    *np.square(corngrph.dark_field)/ref_intnsty)
#    hci.imshow_field(np.log10(contrast_img))
    hci.imshow_field(contrast_img)
    plt.colorbar()
    plt.show()
    contrast = contrast_img.mean()
    return contrast


def cost_ftn(fpm_wid, lyot_wid, corngrph, incdnt_wavfrnt, ref_intnsty
             , planet_star_ratio=1e-10):
    update_fpm(fpm_wid, corngrph)
    update_lyot(lyot_wid, corngrph)
    prop_output = corngrph(incdnt_wavfrnt)
    contrast = mean_contrast(corngrph, prop_output, ref_intnsty)
    peak_factor = ( fpm_thrghpt(prop_output, corngrph)*lyot_thrghpt(corngrph)
                                          *strehl(prop_output, ref_intnsty)
                  )
    print("Contrast = {:E}".format(contrast)) 
    print("Peak factor = {:E}".format(peak_factor)) 
    return (contrast+planet_star_ratio)/peak_factor


def cost_ftn_2(fpm_wid, lyot_wid, corngrph, incdnt_wavfrnt, ref_intnsty
             , planet_star_ratio=1e-10):
    update_fpm(fpm_wid, corngrph)
    update_lyot(lyot_wid, corngrph)
    prop_output = corngrph(incdnt_wavfrnt)
    contrast = mean_contrast(corngrph, prop_output, ref_intnsty)
    raw_cost = (contrast+planet_star_ratio)/( (fpm_thrghpt(prop_output, corngrph)
                                          *lyot_thrghpt(corngrph)
                                          *strehl(prop_output, ref_intnsty)) )
    beta = 1.0
    target = 1e-10
    tolrnc = 3e-10
    cost = -beta*((raw_cost-target)/tolrnc)**2
    return cost


def demo():
    t0 = time.time()   
    coro = ini_corngrph()
    output = coro(hci.Wavefront(coro.pupil_field, coro.wvl))
    print("ini_corngrph & propagation took {} s".format(time.time()-t0))
    plt.subplot(221)
    hci.imshow_field(coro.fpm_field, cmap='inferno')
    plt.colorbar()
    plt.subplot(222)
    hci.imshow_field(output['camera'].intensity, cmap='inferno')
    plt.colorbar()
    t0 = time.time()   
    update_fpm(3.0*coro.f_num*coro.wvl, coro)
    output = coro(hci.Wavefront(coro.pupil_field, coro.wvl))
    print("update_fpm and propagation took {} s".format(time.time()-t0))
    plt.subplot(223)
    hci.imshow_field(coro.fpm_field, cmap='inferno')
    plt.colorbar()
    plt.subplot(224)
    hci.imshow_field(output['camera'].intensity, cmap='inferno')
    plt.colorbar()
    plt.show()


def _test_fpm_thrghpt():
    """Test <fpm_thrghpt> function.  To get valid results, 
    fpm_coef must be > 2*iwa and owa <= 28 (outer extent of focal grid)"""

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
    hci.imshow_field(output0['on_fpm'].intensity, cmap='inferno')
    plt.colorbar()
    plt.subplot(222)
    hci.imshow_field(output0['fpm'].intensity, cmap='inferno')
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


def _test_cost_ftn():
    coro = ini_corngrph()    
    fpm_wid = 5.0*coro.wvl/coro.pupil_diamtr
    lyot_wid = 0.6*coro.pupil_diamtr
    wf = hci.Wavefront(coro.pupil_field, coro.wvl)
    output = coro(wf)
    ref_intnsty = output['on_fpm'].intensity.max()
    result = cost_ftn(fpm_wid, lyot_wid, coro, wf, ref_intnsty)
    print("Cost function value = {}".format(result))


def _try_mystic():
    from mystic.math import almostEqual
    from mystic.monitors import Monitor 
    from mystic.solvers import NelderMeadSimplexSolver  
    from mystic.termination import CandidateRelativeTolerance as Crt  
    # Set up coronagraph     
    print("Setting up coronagraph")
    coro = ini_corngrph(input_args())
    wf = hci.Wavefront(coro.pupil_field, coro.wvl)
    output = coro(wf)
    ref_intnsty = output['on_fpm'].intensity.max()
    # Set optimization parameters
    print("Initializing optimization parameters")
    maxiter = 100 
    print("maxiter = {}".format(maxiter))
    ndim = 2 
    print("ndim = {}".format(ndim))
    x0 = [5.0*coro.wvl/coro.pupil_diamtr, 0.7*coro.pupil_diamtr]
    print("x0 = {}".format(x0))
    lower_bound = [1.0*coro.wvl/coro.pupil_diamtr, 0.5*coro.pupil_diamtr]
    print("lower_bound = {}".format(lower_bound))
    upper_bound = [8*coro.wvl/coro.pupil_diamtr, 1.0*coro.pupil_diamtr]
    print("upper_bound = {}".format(upper_bound))
    esow = Monitor() 
    ssow = Monitor()
    solver = NelderMeadSimplexSolver(ndim)
    solver.SetInitialPoints(x0)
    solver.SetStrictRanges(lower_bound, upper_bound)
    solver.SetEvaluationLimits(generations=maxiter)
    solver.SetEvaluationMonitor(esow)
    solver.SetGenerationMonitor(ssow)
    term = Crt()
    time1 = time.time() 
    print("wall time = {}".format(time1))
    # Set up cost function
    "Setting up cost function"
    def costfunc(x):
        return cost_ftn(x[0], x[1], coro, wf, ref_intnsty)
    # Start solving 
    "Start solving"
    solver.Solve(costfunc, term)
    sol = solver.Solution()
    time_elapsed = time.time() - time1
    fx = solver.bestEnergy
    output_final = coro(wf)
    contrast_final = output_final['camera'].intensity/ref_intnsty
    print("Solution:  {}".format(sol))
    print("Solution {}*wvl/D, {}*D".format(sol[0]/(coro.wvl/coro.pupil_diamtr)
                                           , sol[1]/coro.pupil_diamtr))
    print("f value:  {}".format(fx))
    print("Iterations:  {}".format(solver.generations))
    print("Function evaluations:  {}".format(len(esow.x)))
    print("Time elapsed:  {} s".format(time_elapsed))
    print("Avg. raw contrast:  {}".format(mean_contrast(coro, output_final
                                          , ref_intnsty)))
    hci.imshow_field(np.log10(contrast_final), vmin='-12', vmax='0'
                     , cmap='cubehelix')
    plt.colorbar()
    plt.show()
    return sol, esow, ssow, coro, wf, ref_intnsty

    #assert almostEqual(fx, 0.0, tol=1e-11)


def _test_contrast():
    coro = ini_corngrph(input_args(oversmpl=4, num_airy=128, field_stop=True))   
    wf = hci.Wavefront(coro.pupil_field, coro.wvl)
    output = coro(wf)
    ref_intnsty = output['on_fpm'].intensity.max()
    print("Contrast = {}".format(mean_contrast(coro, output, ref_intnsty)))


if __name__ == '__main__':
#    _try_mystic()
    _test_contrast()
