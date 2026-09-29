""" sleec04.py -- Manipulates the optical elements of coro_sleec.Coronagraph """

import os
import numpy as np
import hcipy_modules as hci
import sleec_optics as sop
from  scipy.optimize import minimize
import matplotlib.pyplot as plt
from   matplotlib.colors import ListedColormap
import time
import fits3
import turbo_colormap as tclr
import mkdonut as mkd
import yaml

def input_args():
    # Coronagraph configuration
    cfg = yaml.full_load(open('sleecParms.yml','r'))

    pars = {}                   # create coro dictionary   
    mthd = cfg['method']        # scipy.optimize.minimize method
    pars['mthd'] = mthd
    wvl = cfg['wavelength']
    pars['wvl'] = float(wvl)
    f_num = cfg['f_number']
    pars['f_num'] = float(f_num)
    flamD = wvl * f_num
    num_acttrs_across = cfg['pupil_dia_actuators']  # Number of actuators across the DM aperture
    acttr_pxl_spacing = cfg['pixels_per_actuator']  # Pixel sampling of the super-influence function
    acttr_spacing = cfg['actuator_pitch']    # Actuator pitch (meters)
    num_pxl_across = 2*(num_acttrs_across*acttr_pxl_spacing + 1)  # Number of pixels across the pupil
    pupil_diamtr = num_acttrs_across*acttr_spacing  # Physical pupil diameter at DM1

    # `pupil_grid`
    pars['pupil_diamtr'] = float(pupil_diamtr)
    pupil_grid = hci.make_pupil_grid(dims=num_pxl_across, diameter=pupil_diamtr)

    # `pupil_field` 
    pupil_field = hci.circular_aperture(diameter=0.98*pupil_diamtr)(pupil_grid)
    pars['pupil_field'] = pupil_field
        
    ### Free parameters ###
    p0 = []
    bnds = []

    ### Focal Plane Mask ###
    pp,b1,b2 = np.multiply(cfg['fpm_diameter'],flamD);      p0.append(pp); bnds.append([b1/pp,b2/pp])
    pp,b1,b2 = cfg['fpm_od'];                               p0.append(pp); bnds.append([b1/pp,b2/pp])
    pp,b1,b2 = np.multiply(cfg['fpm_cosine_apod'],flamD);   p0.append(pp); bnds.append([b1/pp,b2/pp])
    pp,b1,b2 = np.multiply(cfg['dielectric_dia'],flamD);    p0.append(pp); bnds.append([b1/pp,b2/pp])
    pp,b1,b2 = np.multiply(cfg['dielectric_phase'],np.pi);  p0.append(pp); bnds.append([b1/pp,b2/pp])
    pp,b1,b2 = np.multiply(cfg['dielectric_apod'],flamD);   p0.append(pp); bnds.append([b1/pp,b2/pp])
    fpm_occ = hci.complex_circular_occulter(p0[0:7])

    ### Lyot Mask ###
    pp,b1,b2 = np.multiply(cfg['lyot_diameter'],pupil_diamtr);  p0.append(pp); bnds.append([b1/pp,b2/pp])
    lyot_dia = pp
    lyot_transm = hci.circular_aperture(lyot_dia)(pupil_grid)
    pars['lyot_field'] = lyot_transm

    ### DM wtransmitted wavefront ###
    pp,b1,b2 = cfg['dm1_z04'];      p0.append(pp); bnds.append([b1,b2])
    pp,b1,b2 = cfg['dm1_z11'];      p0.append(pp); bnds.append([b1,b2])
    pp,b1,b2 = cfg['dm1_z22'];      p0.append(pp); bnds.append([b1,b2])
    pp,b1,b2 = cfg['dm1_z37'];      p0.append(pp); bnds.append([b1,b2])
    pp,b1,b2 = cfg['dm2_z04'];      p0.append(pp); bnds.append([b1,b2])
    pp,b1,b2 = cfg['dm2_z11'];      p0.append(pp); bnds.append([b1,b2])
    pp,b1,b2 = cfg['dm2_z22'];      p0.append(pp); bnds.append([b1,b2])
    pp,b1,b2 = cfg['dm2_z37'];      p0.append(pp); bnds.append([b1,b2])

    # collect free parameters, bounds, and contraints
    pars['p0'] = p0
    pars['bnds'] = bnds

    ## contraints ###
    cons = []
    e1 = {'type' : 'eq', 'fun' : lambda x, : x[0] - x[3]}      # link parameters as constraints for SLSQP and COBYLA
    e2 = {'type' : 'eq', 'fun' : lambda x, : x[2] - x[5]}      # link parameters as constraints for SLSQP and COBYLA
    cons.append(e1)
    cons.append(e2)

    if mthd == 'COBYLA':                                      # cast bounds as constraints for COBYLA only
        for nn in range(len(bnds)):
            lower, upper = bnds[nn]
            l = {'type' : 'ineq', 'fun' : lambda x, lb=lower, i=nn : x[i] - lb}
            u = {'type' : 'ineq', 'fun' : lambda x, ub=upper, i=nn : ub - x[i]}
            cons.append(l)
            cons.append(u)

    pars['cons'] = cons

    # return to coro configuration
    
    # `fpm_grid` 
    pix_per_airy = 16   # number of pixels per (wvl * f_num) 
    pars['pix_per_airy'] = pix_per_airy
    num_airy = 24     # Outer extent of the computational field in terms of Airy rings
    focal_grid = hci.make_focal_grid(pupil_grid=pupil_grid, q=pix_per_airy, num_airy=num_airy,
        focal_length=f_num*pupil_diamtr, wavelength=wvl)

    # `fpm_field` 
    fpm_transm = fpm_occ(focal_grid) 
    pars['fpm_field'] = fpm_transm

    #  `dm1_field' (in meters)
    pars['dm1_surfc'] = hci.Field(np.zeros(num_pxl_across**2),pupil_grid)        
    
    #  `dm2_field` (in meters)
    pars['dm2_surfc'] = hci.Field(np.zeros(num_pxl_across**2),pupil_grid)

    # `inflnc_ftn'
    pars['inflnc_ftn'] = None

    #  `dm1_vmap`
    pars['dm1_vmap'] = None
    
    #  `dm2_vmap`
    pars['dm2_vmap'] = None

    ### Fixed parameters ###
   
    #  `dm_seprtn`
    pars['dm_seprtn'] = float(1.0)
    
    # 'IWA' & 'OWA'
    iwa_airy = cfg['IWA']
    pars['iwa_airy'] = float(iwa_airy)
    owa_airy = cfg['OWA']
    pars['owa_airy'] = float(owa_airy)

    # 'contrast_floor'
    floor = cfg['contrast_floor']
    pars['contrast_floor'] = float(floor)

    return pars
    
def ini_coro(args=input_args()):
    # Initialize a coronagraph
    coro = sop.Coronagraph()
    coro.wvl = args['wvl']
    coro.pupil_diamtr = args['pupil_diamtr']
    coro.pupil_field = args['pupil_field']
    coro.fpm_field = args['fpm_field']
    coro.lyot_field = args['lyot_field']
    coro.dm1_surfc = args['dm1_surfc']
    coro.dm2_surfc = args['dm2_surfc']
    coro.dm1_vmap = args['dm1_vmap']
    coro.dm2_vmap = args['dm2_vmap']
    coro.inflnc_ftn = args['inflnc_ftn']
    coro.dm_seprtn = args['dm_seprtn']
    coro.f_num = args['f_num']
    coro.iwa_airy = args['iwa_airy']
    coro.owa_airy = args['owa_airy']
    coro.pix_per_airy = args['pix_per_airy']
    coro.mthd = args['mthd']
    coro.p0 = args['p0']
    coro.bnds = args['bnds']
    coro.cons = args['cons']
    coro.contrast_floor = args['contrast_floor']
    coro.load_elements()
    return coro

def update_pupil(coro, y_tilt_airy, x_tilt_airy):
    ytilt = 2.0*np.pi*y_tilt_airy   # tilt in radians 
    xtilt = 2.0*np.pi*x_tilt_airy   # tilt in radians 
    pupil = hci.tilted_circular_aperture(diameter=0.98*coro.pupil_diamtr, ytilt=ytilt, xtilt=xtilt)(coro.pupil_field.grid)
    coro.pupil_field = pupil
    coro.optical_elements[0] = coro.pupil
    return coro
    
def update_dm1(xp, coro):
    zs = np.zeros(37)   # C-index of zernike terms, hence focus z(4) is indexed as z[3]
    zs[3]  = xp[7]   # z(4) focus
    zs[10] = xp[8]   # z(11) spherical
    zs[21] = xp[9]
    zs[36] = xp[10]
    dm1_surfc = hci.zernike_aperture(zs,diameter=0.98*coro.pupil_diamtr)(coro.pupil_field.grid)
    coro.dm1_surfc = dm1_surfc
    coro.optical_elements[1] = coro.dm1
    return coro

def update_dm2(xp, coro):
    zs = np.zeros(37)
    zs[3]  = xp[11]
    zs[10] = xp[12]
    zs[21] = xp[13]
    zs[36] = xp[14]
    dm2_surfc = hci.zernike_aperture(zs, diameter=0.98*coro.pupil_diamtr)(coro.pupil_field.grid)
    coro.dm2_surfc = dm2_surfc
    coro.optical_elements[3] = coro.dm2
    return coro

def update_fpm(xp, coro):
    fpm_args = xp[0:6] 
    fpm_occ = hci.complex_circular_occulter(fpm_args)
    fpm_transm = fpm_occ(coro.fpm_field.grid)
    coro.fpm_field = fpm_transm
    coro.optical_elements[6] = coro.fpm
    return coro

def update_lyot(xp, coro):
    lyot_wid = xp[6]
    lyot_transm = hci.circular_aperture(lyot_wid)
    coro.lyot_field = lyot_transm(coro.pupil_field.grid)
    coro.optical_elements[8] = coro.lyot
    return coro

def sum_dark_annulus(inner, outer, coro, prop_output):
    # used for scoring contrast
    donut = ( hci.circular_aperture(2.0*outer)(coro.fpm_field.grid) 
            - hci.circular_aperture(2.0*inner)(coro.fpm_field.grid) )
    sumSpeckle = np.sum(donut*prop_output['camera'].intensity)
    pixSpeckle = np.sum(donut)
    return pixSpeckle, sumSpeckle

def starPSF(coro, incdnt_wavfrnt): 
    # star PSF: returns characteristics of the star PSF
    update_pupil(coro, y_tilt_airy=0.0, x_tilt_airy=0.0)
    output = coro(incdnt_wavfrnt)
    starMax = output['on_fpm'].intensity.max()
    star = output['on_fpm'].intensity
    ii = np.where(star >= 0.5*starMax)
    mask = np.zeros_like(star)
    mask[ii] = 1
    pixStarPSF = np.sum(mask)
    sumStarPSF = np.sum(mask*star)
#    print ('reference star: %6.3e  %6.3e  %5.0f' % (starMax, sumStarPSF, pixStarPSF))
    return starMax, sumStarPSF, pixStarPSF

def planetPSF(coro, incdnt_wavfrnt, starMax, sumStarPSF, tilt_airy):
    # planet PSF returns number of pixels within the FWHM  and pixel intensities summed 
    #     over the FWHM, for the offset star in the coro (camera) focal plane
    update_pupil(coro, y_tilt_airy=tilt_airy, x_tilt_airy=0.0)
    output_tilted = coro(incdnt_wavfrnt)
    PSF = output_tilted['camera'].intensity
    mask = np.zeros_like(PSF)
    ii = np.where(PSF >= 0.5*PSF.max())
    mask[ii] = 1
    sumPSF = np.sum(mask*PSF)
    pixPSF = np.sum(mask)
    hwhmPSF = np.sqrt(pixPSF/np.pi) / coro.pix_per_airy  # define annular scoring footprint
    r1 = (tilt_airy - hwhmPSF)*coro.f_num*coro.wvl
    r2 = (tilt_airy + hwhmPSF)*coro.f_num*coro.wvl
    update_pupil(coro, y_tilt_airy=0.0, x_tilt_airy=0.0)
    output_untilted = coro(incdnt_wavfrnt)
    pixDonut, sumDonut = sum_dark_annulus(r1, r2, coro, output_untilted)
    makefits = False
    if makefits:
        imgdim = int(np.sqrt(PSF.shape))
        img = np.reshape(PSF,[imgdim,imgdim])
        fits3.wfits("planet_image.fits",img)
#        fits3.wfits("donut_image.fits",img+dmask*img)
    return  pixPSF, sumPSF, pixDonut, sumDonut

def cost_ftn(x, args):
    coro = args[0]
    wf = args[1]
    p0 = args[2]
    starMax = args[3]
    sumStarPSF = args[4]
    contrast_floor = args[5]
    xp = x * p0
    update_fpm(xp, coro)
    update_lyot(xp, coro)
    update_dm1(xp, coro)
    update_dm2(xp, coro)
#    print( '%8.3f  %6.3f  %8.3f  %8.3f  %6.3f  %6.3f' % (xp[0]*1e6, xp[1], xp[2]*1e6, xp[3]*1e6, xp[4], xp[5]*1e6))
    print( '%8.3f  %6.3f  %8.3f  %8.3f  %6.3f  %6.3f  %6.3f  %6.3e  %6.3e  %6.3e  %6.3e  %6.3e  %6.3e  %6.3e  %6.3e' % 
            (xp[0]*1e6, xp[1], xp[2]*1e6, xp[3]*1e6, xp[4], xp[5]*1e6, xp[6], xp[7], xp[8], xp[9], xp[10], xp[11], xp[12], xp[13], xp[14]))
    tilt_airy = 3.5                                                           # thruput is estimated for a planet PSF at this offset from star
    pixPlntPSF, sumPlntPSF, pixDonut, sumDonut = planetPSF(coro, wf, starMax, sumStarPSF, tilt_airy)
    thruput = sumPlntPSF / sumStarPSF
    mean_speckle = sumDonut / pixDonut                           # mean speckle intensity per pixel at the annular field centered on the planet PSF footprint 
    mean_contrast = mean_speckle / (thruput*starMax)             # standard definition of contrast, pixel sampling of dark field and star PSF
#    raw_cost = (mean_contrast * pixPlntPSF) * (starMax/sumStarPSF) / thruput   # (dark field sum / star sum) over core PSF at coro focal plane / thruput 
    raw_cost = (mean_contrast + contrast_floor) * pixPlntPSF / thruput                      # same, but factor out constant terms
    print ('pixPlntPSF=%8.2f sumPlntPSF=%6.3e  pixDonut=%8.2f  sumDonut=%6.3e    cost=%7.3e   mean_contrast=%7.3e   thruput=%6.3f' %
             (pixPlntPSF, sumPlntPSF, pixDonut, sumDonut, raw_cost, mean_contrast, thruput))
    return raw_cost

def cost_final(xp, coro, wf, starMax, sumStarPSF, contrast_floor):    # same as cost_fcn() for final report
    update_fpm(xp, coro)
    update_lyot(xp, coro)
    update_dm1(xp, coro)
    update_dm2(xp, coro)
    tilt_airy = 3.5                                                           # thruput is estimated for a planet PSF at this offset from star
    pixPlntPSF, sumPlntPSF, pixDonut, sumDonut = planetPSF(coro, wf, starMax, sumStarPSF, tilt_airy)
    thruput = sumPlntPSF / sumStarPSF
    mean_speckle = sumDonut/ pixDonut                                         # mean speckle intensity per pixel at the coro focal plane 
    mean_contrast = mean_speckle / (thruput*starMax)                          # standard definition of contrast
    raw_cost = (mean_contrast + contrast_floor) * pixPlntPSF / thruput        # factor out constant terms
    return raw_cost, mean_contrast, thruput

def _try_optimize():

    # Set up coronagraph     
    coro = ini_coro(input_args())
    wf = hci.Wavefront(coro.pupil_field, coro.wvl)
    output = coro(wf)
    starMax, sumStarPSF, pixStarPSF = starPSF(coro, wf)  # this does not change -- compute once only.
    
    # Set optimization parameters
    maxitr = 1000
    flamD = coro.f_num * coro.wvl
    pupdia = coro.pupil_diamtr
    contrast_floor = coro.contrast_floor

    p0 = coro.p0
    x0 = np.ones(len(p0))           # preconditioned parameters all start out with value = 1.0
    bnds = coro.bnds               # bounds are ratios to the intial parameters
    cons = coro.cons               # constraints are all ratios to the initial parameters

    # Solve 
    scimthd = coro.mthd
    arg = []
    arg.append(coro)
    arg.append(wf)
    arg.append(p0)
    arg.append(starMax)
    arg.append(sumStarPSF)
    arg.append(contrast_floor)
    time1 = time.time() 

    # Nelder-Mead
    if scimthd == 'Nelder':
        result = minimize(cost_ftn, x0, args=(arg), method='Nelder-Mead', tol=1e-5, options={'maxiter' : maxitr})
    
    # L-BFGS-B = Limited-memory Broyden-Fletcher-Goldfarb-Shanno algorithm with box constraints
    if scimthd == 'L-BFGS-B':
        result = minimize(cost_ftn, x0, args=(arg), method='L-BFGS-B', bounds=bnds, options={'maxiter' : maxitr, 'eps' : 0.001, 'ftol' : 1e-10})
#        result = fmin_l_bfgs_b(cost_ftn, x0, args=(arg), bounds=bnds, approx_grad=True, epsilon=0.01, maxiter=maxitr, factr=10)

    # TNC = truncated Newton method
    if scimthd == 'TNC':
        result = minimize(cost_ftn, x0, args=(arg), method='TNC', tol=1e-11, bounds=bnds, options={'maxiter' : maxitr, 'ftol' : 1e-5, 'eps' : 0.00001})
#        result = minimize(cost_ftn, x0, args=(arg), method='TNC', tol=1e-10, options={'maxiter' : maxitr, 'eps' : 0.001})

    # COBYLA = Constained Optimization BY Constrained Approximation 
    if scimthd == 'COBYLA':
        result = minimize(cost_ftn, x0, args=(arg), method='COBYLA', constraints=cons, options={'rhobeg' : 0.01, 'maxiter' : maxitr, 'tol' : 0.00001})
#        result = minimize(cost_ftn, x0, args=(arg), method='COBYLA', options={'rhobeg' : 0.0001, 'maxiter' : maxitr, 'tol' : 0.0001})

    # SLSQP = Sequential Least SQuares Programming
    if scimthd == 'SLSQP':
        result = minimize(cost_ftn, x0, args=(arg), method='SLSQP', constraints=cons, bounds=bnds, options={'maxiter' : maxitr, 'ftol' : 1e-11, 'eps' : 0.0001})


    time_elapsed = time.time() - time1
    print ('Final Cost: %7.3e' % result.fun)
    if scimthd == 'COBYLA':
        print ('Function evaluations: %4d' % (result.nfev))
    else:
        print ('Iterations: %4d, Function evaluations: %4d' % (result.nit, result.nfev))
    print ('Final Status:', result.status)
    print ('Message:', result.message)
 
    output_final = coro(wf)
    coro_final = output_final['camera'].intensity/starMax
    flamD = coro.f_num*coro.wvl
    xp = result.x * p0
    eminutes = time_elapsed/60.0
    print("Solution: FPM diameter =%6.3f*f_num*wvl, FPM OD =%6.3f, FPM apodization =%6.3f*f_num*wvl"
            % (xp[0]/flamD, xp[1], xp[2]/flamD) )
    print("Solution: Dielectric dia =%6.3f*f_num*wvl, Phase shift =%6.3f*pi, Dielectric apod =%6.3f*f_num*wvl"
            % (xp[3]/flamD, xp[4]/np.pi, xp[5]/flamD) )
    print("Solution: Lyot mask diameter =%6.3f*D" % (xp[6]/(coro.pupil_diamtr)))
    print('Time elapsed: %5.1f minutes' % (eminutes))

    print()
    print(result)

    cost, contrast, thruput = cost_final(xp, coro, wf, starMax, sumStarPSF, contrast_floor)    # same as cost_fcn() for final report
    print('result: ', cost, contrast, thruput)

    # archive the results
    resdir = '../results/'
    timestamp = time.strftime('%Y%m%d_%H%M',time.localtime())

    kwd = {}                                     # save results as FITS keywords
    for k in range(len(xp)):
        knme = 'xinit%02d' % k
        kwd[knme] = p0[k]
        knme = 'xfinal%02d' % k
        kwd[knme] = xp[k]
    kwd['floor'] = float(contrast_floor)
    kwd['cost'] = float(cost)
    kwd['contrast'] = float(contrast)
    kwd['thruput'] = float(thruput)
    if scimthd != 'COBYLA':
        kwd['iters'] = result.nit
    else:
        kwd['iters'] = 0.0
    kwd['evals'] = result.nfev
    kwd['method'] = scimthd
    kwd['minutes'] = eminutes

    fpm_transm = coro.fpm_field                  # archive complex focal plane mask
    fpmdim = int(np.sqrt(fpm_transm.shape))
    mask = np.reshape(fpm_transm,[fpmdim,fpmdim])
    amask = np.abs(mask)     # amplitude transmission
    aamask = amask*amask     # intensity transmittance
    fname = resdir+'FPM_trans_'+timestamp+'.fits'
    fits3.wfits(fname, aamask, keywordvals=kwd)
    pmask = np.angle(mask)
    fname = resdir+'FPM_phase_'+timestamp+'.fits'
    fits3.wfits(fname, pmask, keywordvals=kwd)

    imgdim = int(np.sqrt(coro_final.shape))       # archive coronagraph image
    img = np.reshape(coro_final,[imgdim,imgdim])
    cimg = imgdim//2 - 0.5
    iwa = coro.iwa_airy * coro.pix_per_airy
    owa = coro.owa_airy * coro.pix_per_airy
    c0 = mkd.mkdonut(img,[cimg,cimg], iwa-2, iwa-1, iwa+1, iwa+2)
    c1 = mkd.mkdonut(img,[cimg,cimg], owa-2, owa-1, owa+1, owa+2)
    cc0 = c0.reshape([imgdim*imgdim,])
    cc1 = c1.reshape([imgdim*imgdim,])
    fname = resdir+'coro_image_'+timestamp+'.fits'
    fits3.wfits(fname, img+c0+c1, keywordvals=kwd)

    turbo = ListedColormap(tclr.turbo_colormap_data)
    hci.imshow_field(np.log10(coro_final+cc0+cc1), vmin='-12', vmax='-4', cmap=turbo )
    plt.colorbar()
    plt.show()

    return
if __name__ == '__main__':
    _try_optimize()



