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
    coefs = design_parameters[:(len(input_pars['n2d']))]
    fpm_dia = design_parameters[len(input_pars['n2d'])]
    lyot_dia = design_parameters[len(input_pars['n2d'])+1]
    if (
        (input_pars['optmz_pupil_coefs'][1] <= coefs).all() 
        and (coefs <= input_pars['optmz_pupil_coefs'][2]).all() 
        and (input_pars['optmz_fpm_dia'][1] <= fpm_dia 
            <= input_pars['optmz_fpm_dia'][2])
        and (input_pars['optmz_lyot_dia'][1] <= lyot_dia 
            <= input_pars['optmz_lyot_dia'][2])
       ):
        print("Parameter priors passed")
        return 0.0
    else:
        print("Parameter priors FAILED")
        return -np.inf
 

def ln_probability(design_parameters, coronagraph, input_pars, non_negative
                   , reconstruct=False):
    """
    Return log probability for optimization using EMCEE.  
    
    """
    print("\n\n************\ndesign_parameters:  {}".format(design_parameters))
    coefs = design_parameters[:(len(input_pars['n2d']))]
    coronagraph = sla.pupil_slepian_r(coefs=coefs
                        , corngrph=coronagraph
                        , file_exists=True, filename=input_pars['filename']
                        , non_negative=non_negative
                        , n2d=None, L=None, method=None
                        , orthgnlty=None, num_dense_pts=None)
    coronagraph = sla.fpm_dot(fpm_dia=design_parameters[len(input_pars['n2d'])]
                , corngrph=coronagraph
                , num_airy=input_pars['num_airy'], field_stop=True)
    coronagraph = sla.lyot_circ_clear(
                        lyot_dia=design_parameters[len(input_pars['n2d'])+1]
                        , corngrph=coronagraph)
    if reconstruct:
        wavefronts = mrt[2][3]
        mrt = slm.merit(coronagraph=coronagraph, inner_radius=input_pars['IWA']
                          , outer_radius=input_pars['OWA'], amplitude1=1.0
                          , area_penalty=input_pars['area_penalty']
                          , off_axis_locatn=off_axis
                          , azimuthal_locatn=2*np.pi/7, amplitude2=1.0
                          , flux_ratio=float(input_pars['flux_ratio'])
                          , output_wavefronts=True)
        whole_throughput, throughputs, planet_core_area = mrt[1]
        whole_contrast, contrasts, area_penalty, wavefronts = mrt[2]
        ln_prob = np.log(mrt[0])
        return (ln_prob, whole_throughput, planet_core_area, whole_contrast
                , area_penalty
                , coronagraph, wavefronts, throughputs, contrasts)
    ln_prr = ln_prior(design_parameters, input_pars)
    if ((coronagraph.pupil_field>=0.0).all() 
        and (coronagraph.pupil_field <= 1.0).all()):  
#    if ((coronagraph.pupil_field <= 1.0).all()):  
        print("Pupil field within limits")
        ln_pupil_pr = 0.0
    else:
        print("Pupil field EXCEEDS LIMITS")
        ln_pupil_pr = -np.inf
    if (np.isfinite(ln_prr) and np.isfinite(ln_pupil_pr)):
        off_axis = (input_pars['IWA'] + input_pars['OWA'])/2.0
        mrt = slm.merit(coronagraph=coronagraph, inner_radius=input_pars['IWA']
                          , outer_radius=input_pars['OWA'], amplitude1=1.0
                          , area_penalty=input_pars['area_penalty']
                          , off_axis_locatn=off_axis
                          , azimuthal_locatn=2*np.pi/7, amplitude2=1.0
                          , flux_ratio=float(input_pars['flux_ratio'])
                          , output_wavefronts=False)
        ln_prob = ln_prr + np.log(mrt[0]) + ln_pupil_pr
        whole_throughput, throughputs, planet_core_area = mrt[1]
        throughputs = throughputs.tolist
        whole_contrast, contrasts, area_penalty = mrt[2]
        contrasts = contrasts.tolist
#    if not (np.isfinite(ln_prr) and np.isfinite(ln_pupil_pr)):
    else:
        ln_prob = -np.inf
        whole_throughput, throughputs, planet_core_area = [0.0, 0.0, 0.0]
        whole_contrast, contrasts, area_penalty = [0.0, 0.0, 0.0]
#    else:
#        ln_prob = ln_prr + np.log(mrt[0]) + ln_pupil_pr
    print("*** ln probability = {} ***".format(ln_prob))
    return (ln_prob, whole_throughput, planet_core_area
            , whole_contrast, area_penalty)
    

def run_emcee(non_negative, spread_factor=0.001, parallel=False):
    with open('sleec_pars.yml', 'r') as f:
        input_pars = yml.full_load(f)
    print("Loaded input parameters from user configuration file")
    save = input("Save EMCEE run to a file?  [y/n]:  ")
    if save == 'y':
        fn_prefix = input("Enter filename prefix for saving EMCEE run:  ")
        path = os.path.join('..', '..', 'sleec_out')
        file_path_hdf = os.path.join(path, fn_prefix+'.hdf')
        file_path_slp = os.path.join(path, fn_prefix+'.pkl')
        file_path_yml = os.path.join(path, fn_prefix+'.yml')
        backend = emc.backends.HDFBackend(file_path_hdf)
        backend.reset(input_pars['nwalkers'], input_pars['ndim'])
        stl.copyfile(input_pars['filename'], file_path_slp)    
        stl.copyfile('sleec_pars.yml', file_path_yml)    
    else:
        backend = None
    coronagraph = sli.import_pars()
    output0 = coronagraph(hci.Wavefront(coronagraph.pupil_field
                                        , coronagraph.wvl))
    print("Initialized coronagraph")
    sla.pupil_slepian_r(coefs=np.ones(len(input_pars['n2d']))
                        , corngrph=coronagraph, file_exists=False
                        , filename=input_pars['filename']
                        , non_negative=True
                        , n2d=input_pars['n2d'], L=170
                        , method='SE', orthgnlty=0, num_dense_pts=4096)
    print("Slepian file created")
    ini_positns = np.empty((input_pars['nwalkers'], input_pars['ndim']))

    for i in range(len(input_pars['n2d'])):
        par_range = spread_factor*(input_pars['optmz_pupil_coefs'][2]
                          - input_pars['optmz_pupil_coefs'][1])
        ini_positns[:, i] = ( 
                             input_pars['optmz_pupil_coefs'][0] 
                             + par_range
                               *np.random.randn(input_pars['nwalkers'])                             )
    ini_positns[:, len(input_pars['n2d'])] = ( 
        input_pars['optmz_fpm_dia'][0] 
        + spread_factor*(input_pars['optmz_fpm_dia'][2]
                         -input_pars['optmz_fpm_dia'][1])
          *np.random.randn(input_pars['nwalkers'])
                                   ) 
    ini_positns[:, len(input_pars['n2d'])+1] = ( 
        input_pars['optmz_lyot_dia'][0] 
        + spread_factor*(input_pars['optmz_lyot_dia'][2]
        -input_pars['optmz_lyot_dia'][1])
          *np.random.randn(input_pars['nwalkers'])
                                   ) 
#    for q in range(input_pars['nwalkers']):
#        for i in range(len(input_pars['n2d'])):
#            maximum = input_pars['optmz_pupil_coefs'][2]
#            minimum = input_pars['optmz_pupil_coefs'][1]
#            ini_positns[q, i] = np.random.uniform(minimum, maximum)
#        ini_positns[q, len(input_pars['n2d'])] = ( 
#            np.random.uniform(input_pars['optmz_fpm_dia'][1]
#                              , input_pars['optmz_fpm_dia'][2])
#                                          )
#        ini_positns[q, len(input_pars['n2d'])+1] = ( 
#            np.random.uniform(input_pars['optmz_lyot_dia'][1]
#                              , input_pars['optmz_lyot_dia'][2])
#                                            )
    print("Initial parameter values (walker positions) set")
    print(ini_positns)  # debug
    dtype = [("whole_throughput", float) 
             , ("planet_core_area", float), ("whole_contrast", float)
             , ("area_penalty", float)]

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
                                         , non_negative
                                                 )
                                         , blobs_dtype=dtype)
         print("Ensemble Sampler instantiated.  Starting optimization run now.")
         sampler.run_mcmc(ini_positns, input_pars['nsteps'], progress=True
                          , skip_initial_state_check=False)


if __name__ == '__main__':  
    run_emcee(non_negative=True, spread_factor=0.001, parallel=True)






