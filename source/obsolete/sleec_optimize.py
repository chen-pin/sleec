import numpy as np


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
    """
    Return average contrast in a zone defined by <corngrph.dark_field>, whose 
    elements are averaging weights.  
    """
    contrast_img = (prop_output['camera'].intensity
                    *corngrph.dark_field/ref_intnsty)
    contrast = contrast_img.sum()/corngrph.dark_field.sum()
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
    hci.imshow_field(output0['on_fpm'].intensity, cmap=turbo)
    plt.colorbar()
    plt.subplot(222)
    hci.imshow_field(output0['fpm'].intensity, cmap=turbo)
    plt.colorbar()
    plt.subplot(223)
    hci.imshow_field(corngrph.dark_field, cmap=turbo)
    plt.colorbar()
    plt.subplot(224)
    hci.imshow_field(corngrph.fpm.mask_transm, cmap=turbo)
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
    hci.imshow_field(np.absolute(corngrph.pupil.mask_transm)**2, cmap=turbo)
    plt.colorbar()
    plt.subplot(122)
    hci.imshow_field(np.absolute(corngrph.lyot.mask_transm)**2, cmap=turbo)
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
    
    
