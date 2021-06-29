"""<sleec_run>:  Module for running optimizations."""
import os
import shutil as stl
import pickle as pkl
from multiprocessing import Pool
import yaml as yml
import numpy as np
import matplotlib.pyplot as plt
import emcee as emc
import hcipy as hci
import sleec_input as sli
import sleec_alter as sla
import sleec_merit as slm


def ln_prior(design_parameters, input_pars):
    """
    Retrun log probability of design parameters
    
    """
    count = 0
    for par in design_parameters:
        if ( par < input_pars['bounds_pars'][count][1] 
             or par > input_pars['bounds_pars'][count][2] ):
            print("Parameter {} out of bound".format(count)) 
            return -np.inf
        count += 1
    return 0.0
 

def ln_probability(design_parameters, coronagraph, input_pars, non_negative
                   , reconstruct=False):
    """
    Return log probability for optimization using EMCEE.  
    
    """
    print("\n\n************\ndesign_parameters:  {}".format(design_parameters))
    pupil_pars = design_parameters[:input_pars['num_pupil_pars']]
    if input_pars['pupil_type'] == 'slepian_n2d':
        coronagraph = sla.pupil_slepian_r(coefs=pupil_pars
                          , corngrph=coronagraph
                          , file_exists=True
                          , filename=input_pars['pupil_filename']
                          , non_negative=non_negative
                          , n2d=None, L=None, method=None
                          , orthgnlty=None, num_dense_pts=None)
    ndx0 = input_pars['num_pupil_pars']
    ndx1 = ndx0 + input_pars['num_dm1_pars']
    dm1_pars = design_parameters[ndx0 : ndx1]
    if input_pars['dm1_type'] == 'disk_harmonic':
        coronagraph = sla.dm_disk_harmonic(coefs=dm1_pars
                          , corngrph=coronagraph
                          , dm=input_pars['dm_dm1']
                          , file_exists=True
                          , filename=input_pars['dm1_filename']
                          , m=None
                          , n_min=None
                          , n_max=None
                          , bc=None)
    elif input_pars['dm1_type'] == 'slepian_fourier':
        coronagraph = sla.dm1_slepian_2d(coefs=dm1_pars
                          , corngrph=coronagraph
                          , file_exists=True
                          , filename=input_pars['dm1_filename']
                          , m_max=None
                          , n2d=None
                          , k=None
                          , L=None
                          , method=None
                          , orthgnlty=None, num_dense_pts=None)
    ndx0 = ndx1
    ndx1 = ndx0 + input_pars['num_dm2_pars']
    dm2_pars = design_parameters[ndx0 : ndx1]
    if input_pars['dm2_type'] == 'disk_harmonic':
        coronagraph = sla.dm_disk_harmonic(coefs=dm2_pars
                          , corngrph=coronagraph
                          , dm=input_pars['dm_dm2']
                          , file_exists=True
                          , filename=input_pars['dm2_filename']
                          , m=None
                          , n_min=None
                          , n_max=None
                          , bc=None)
    elif input_pars['dm2_type'] == 'slepian_fourier':
        coronagraph = sla.dm2_slepian_2d(coefs=dm2_pars
                          , corngrph=coronagraph
                          , file_exists=True
                          , filename=input_pars['dm2_filename']
                          , m_max=None
                          , n2d=None
                          , k=None
                          , L=None
                          , method=None
                          , orthgnlty=None, num_dense_pts=None, oppose_dm1=True)
    ndx0 = ndx1
    ndx1 = ndx0 + input_pars['num_fpm_pars']
    fpm_pars = design_parameters[ndx0 : ndx1]
    if input_pars['fpm_type'] == 'dot':
        coronagraph = sla.fpm_dot(
                          fpm_dia=fpm_pars
                          , corngrph=coronagraph
                          , num_airy=input_pars['num_airy'])
    elif input_pars['fpm_type'] == 'super_gaussian':
        coronagraph = sla.fpm_supergauss(A=fpm_pars[0], sigma=fpm_pars[1]
                                         , N=fpm_pars[2], corngrph=coronagraph
                                         , num_airy=input_pars['num_airy'])
    ndx0 = ndx1
    ndx1 = ndx0 + input_pars['num_lyot_pars']
    lyot_pars = design_parameters[ndx0 : ndx1]
    if input_pars['lyot_type'] == 'circular_clear':
        coronagraph = sla.lyot_circ_clear(
                          lyot_dia=lyot_pars
                          , corngrph=coronagraph)
    band_centers = np.arange(input_pars['IWA']+0.5, input_pars['OWA']+0.5)
    if reconstruct:
        ( merit_val, tau, throughputs, band_contrasts
             , encircle_psf, wavefronts, Cr_p, Cr_b_star, Cr_b_zodi
             , Cr_b_exozodi, Cr_b ) = slm.merit(
                        coronagraph=coronagraph
                        , weights=input_pars['weights']
                        , inner_radius=input_pars['IWA']
                        , outer_radius=input_pars['OWA']
                        , encircle_radius=input_pars['encircle_radius']
                        , snr=input_pars['snr']
                        , F_0=input_pars['F_0']
                        , m_v=input_pars['m_v']
                        , delta_mag=input_pars['delta_mag']
                        , D=input_pars['D']
                        , wvl=input_pars['wvl']
                        , wvl_bw=input_pars['wvl_bw']
                        , z=input_pars['z']
                        , x=input_pars['x']
                        , n=input_pars['n']    )
        ln_prob = np.log(merit_val)
        whole_throughput = np.sum(throughputs
                                *band_centers)/np.sum(band_centers)
        whole_contrast = np.sum(band_contrasts
                                *band_centers)/np.sum(band_centers)
        return ( coronagraph, merit_val, tau, whole_throughput, throughputs
             , whole_contrast, band_contrasts
             , encircle_psf, wavefronts, Cr_p, Cr_b_star, Cr_b_zodi
             , Cr_b_exozodi, Cr_b ) 
    ln_prr = ln_prior(design_parameters, input_pars)
    if ((coronagraph.pupil_field>=0.0).all() 
        and (coronagraph.pupil_field <= 1.0).all()):  
        print("Pupil field within limits")
        ln_pupil_pr = 0.0
    else:
        print("Pupil field EXCEEDS LIMITS")
        ln_pupil_pr = -np.inf
    if (np.isfinite(ln_prr) and np.isfinite(ln_pupil_pr)):
        ( merit_val, tau, throughputs, band_contrasts
             , encircle_psf, wavefronts, Cr_p, Cr_b_star, Cr_b_zodi
             , Cr_b_exozodi, Cr_b ) = slm.merit(
                        coronagraph=coronagraph
                        , weights=input_pars['weights']
                        , inner_radius=input_pars['IWA']
                        , outer_radius=input_pars['OWA']
                        , encircle_radius=input_pars['encircle_radius']
                        , snr=input_pars['snr']
                        , F_0=input_pars['F_0']
                        , m_v=input_pars['m_v']
                        , delta_mag=input_pars['delta_mag']
                        , D=input_pars['D']
                        , wvl=input_pars['wvl']
                        , wvl_bw=input_pars['wvl_bw']
                        , z=input_pars['z']
                        , x=input_pars['x']
                        , n=input_pars['n']    )
        ln_prob = ln_prr + np.log(merit_val) + ln_pupil_pr
        whole_throughput = np.sum(throughputs
                                *band_centers)/np.sum(band_centers)
        whole_contrast = np.sum(band_contrasts
                                *band_centers)/np.sum(band_centers)
    else:
        ln_prob = -np.inf
        whole_throughput = 0.0
        whole_contrast = 0.0
        encircle_psf = 0.0
    print("*** ln probability = {} ***".format(ln_prob))
    print("avg contrast = {}".format(whole_contrast))
    print("avg throughput = {}".format(whole_throughput))
    return ln_prob, whole_throughput, whole_contrast, encircle_psf
    

def run_emcee(spread_factor=0.001, parallel=False):
    coronagraph, input_pars = sli.import_pars()
    print("Loaded input parameters from user configuration file")
    save = input("Save EMCEE run to a file?  (y/[n]):  ")
    if save == 'y':
        fn_prefix = input("Enter filename prefix for saving EMCEE run:  ")
        path = os.path.join('..', '..', 'sleec_out')
        file_path_hdf = os.path.join(path, fn_prefix+'.hdf')
        file_path_pupil_slp = os.path.join(path, fn_prefix+'_pupil.pkl')
        file_path_dm1_slp = os.path.join(path, fn_prefix+'_dm1.pkl')
        file_path_dm2_slp = os.path.join(path, fn_prefix+'_dm2.pkl')
        file_path_yml = os.path.join(path, fn_prefix+'.yml')
        backend = emc.backends.HDFBackend(file_path_hdf)
        backend.reset(input_pars['nwalkers'], input_pars['ndim'])
        if input_pars['pupil_filename'] != None:
            stl.copyfile(input_pars['pupil_filename'], file_path_pupil_slp)    
        if input_pars['dm1_filename'] != None:
            stl.copyfile(input_pars['dm1_filename'], file_path_dm1_slp)    
        if input_pars['dm2_filename'] != None:
            stl.copyfile(input_pars['dm2_filename'], file_path_dm2_slp)    
        stl.copyfile('sleec_pars.yml', file_path_yml)    
    else:
        backend = None
    non_negative = input_pars['non_negative']
    output0 = coronagraph(hci.Wavefront(coronagraph.pupil_field
                                        , coronagraph.wvl))
    print("Initialized coronagraph")
    if input_pars['pupil_type'] == 'slepian_n2d':
        sla.pupil_slepian_r(coefs=np.ones(input_pars['num_pupil_pars'])
                            , corngrph=coronagraph, file_exists=False
                            , filename=input_pars['pupil_filename']
                            , non_negative=non_negative
                            , n2d=input_pars['n2d_pupil'], L=170
                            , method='SE', orthgnlty=0, num_dense_pts=4096)
        print("Pupil Slepian file created")
    if input_pars['dm1_type'] == 'slepian_fourier':
        sla.dm1_slepian_2d(coefs=np.ones(input_pars['num_dm1_pars'])
                            , corngrph=coronagraph
                            , file_exists=False
                            , filename=input_pars['dm1_filename']
                            , m_max=input_pars['m_dm1']
                            , n2d=input_pars['n2d_dm1']
                            , k=input_pars['k_dm1']
                            , L=170
                            , method='SE'
                            , orthgnlty=0
                            , num_dense_pts=4096)
        print("DM1 Slepian-Fourier file created")
    elif input_pars['dm1_type'] == 'disk_harmonic':
        sla.dm_disk_harmonic(coefs=np.ones(input_pars['num_dm1_pars'])
                            , corngrph=coronagraph
                            , dm='dm1'
                            , file_exists=False
                            , filename=input_pars['dm1_filename']
                            , m=input_pars['m_dm1']
                            , n_min=input_pars['n_min_dm1']
                            , n_max=input_pars['n_max_dm1']
                            , bc=input_pars['bc_dm1'])
        print("DM1 disk-harmoic file created")
    if input_pars['dm2_type'] == 'slepian_fourier':
        sla.dm2_slepian_2d(coefs=np.ones(input_pars['num_dm2_pars'])
                            , corngrph=coronagraph
                            , file_exists=False
                            , filename=input_pars['dm2_filename']
                            , m_max=input_pars['m_dm2']
                            , n2d=input_pars['n2d_dm2']
                            , k=input_pars['k_dm2']
                            , L=170
                            , method='SE'
                            , orthgnlty=0
                            , num_dense_pts=4096)
        print("DM2 Slepian-Fourier file created")
    elif input_pars['dm2_type'] == 'disk_harmonic':
        sla.dm_disk_harmonic(coefs=np.ones(input_pars['num_dm2_pars'])
                            , corngrph=coronagraph
                            , dm='dm2'
                            , file_exists=False
                            , filename=input_pars['dm2_filename']
                            , m=input_pars['m_dm2']
                            , n_min=input_pars['n_min_dm2']
                            , n_max=input_pars['n_max_dm2']
                            , bc=input_pars['bc_dm2'])
        print("DM2 disk-harmoic file created")
    ini_positns = np.empty((input_pars['nwalkers'], input_pars['ndim']))
    i = 0
    for row in input_pars['bounds_pars']:
        par_range = spread_factor*(row[2]-row[1])
        ini_positns[:, i] = (row[0] 
                             + par_range*np.random.randn(input_pars['nwalkers'])
                            )
        i += 1
    print("Initial parameter values (walker positions) set")
    print(ini_positns)  # debug
    dtype = [("whole_throughput", float) 
             , ("whole_contrast", float)
             , ("encircle_psf", float)]

    if parallel:
        os.environ["OMP_NUM_THREADS"] = "1"  # Disable Anaconda's 
            # parallelization
        with Pool() as pool:
            sampler = emc.EnsembleSampler(
                          input_pars['nwalkers'], input_pars['ndim']
                          , ln_probability, backend=backend, pool=pool
                          , args=(coronagraph, input_pars, non_negative)
                          , blobs_dtype=dtype
                                         )
            print("Ensemble Sampler instantiated.")  
            print("Starting optimization run now.")
            sampler.run_mcmc(ini_positns, input_pars['nsteps'], progress=True)
    else:
         sampler = emc.EnsembleSampler(input_pars['nwalkers']
                                         , input_pars['ndim']
                                         , ln_probability, backend=backend 
                                         , args=(coronagraph, input_pars
                                                 , non_negative)
                                         , blobs_dtype=dtype)
         print("Ensemble Sampler instantiated.  Starting optimization run now.")
         sampler.run_mcmc(ini_positns, input_pars['nsteps'], progress=True
                          , skip_initial_state_check=False)


if __name__ == '__main__':  
    run_emcee(spread_factor=0.001, parallel=True)






