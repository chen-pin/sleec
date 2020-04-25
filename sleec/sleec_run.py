"""<sleec_run>:  Module for running optimizations."""
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
    coefs = np.array(design_parameters[0])
    fpm_dia = design_parameters[1]
    lyot_dia = design_parameters[2]
    if (
        input_pars['optmz_pupil_coefs'][1] <= coefs.all() 
            <= input_pars['optmz_pupil_coefs'][2] 
        and input_pars['optmz_fpm_dia'][1] <= fpm_dia 
            <= input_pars['optmz_fpm_dia'][2]
        and input_pars['optmz_lyot_dia'][1] <= lyot_dia 
            <= input_pars['optmz_lyot_dia'][2]):
        return 0.0
    return -np.inf
 

def ln_probability(design_parameters, coronagraph, input_pars):
    """
    Return log probability for optimization using EMCEE.  
    
    """
    print("design_parameters:  {}".format(design_parameters))
    coefs = design_parameters[:(input_pars['k'])]
    sla.pupil_slepian_r(coefs=coefs
                        , corngrph=coronagraph
                        , file_exists=True, filename=input_pars['filename']
                        , m=None, n2d=None, k=None, L=None, method=None
                        , orthgnlty=None, num_dense_pts=None)
    sla.fpm_dot(fpm_dia=design_parameters[input_pars['k']], corngrph=coronagraph
                , num_airy=input_pars['num_airy'], field_stop=True)
    sla.lyot_circ_clear(lyot_dia=design_parameters[input_pars['k']+1]
                        , corngrph=coronagraph)
    off_axis = (input_pars['IWA'] + input_pars['OWA'])/2.0
    mrt = slm.merit(coronagraph=coronagraph, inner_radius=input_pars['IWA']
                      , outer_radius=input_pars['OWA'], amplitude1=1.0
                      , area_penalty=input_pars['area_penalty']
                      , off_axis_locatn=off_axis
                      , azimuthal_locatn=2*np.pi/7, amplitude2=1.0
                      , flux_ratio=float(input_pars['flux_ratio']))[0]
    # Impose limits on the pupil field and FPM field
    if (np.real(coronagraph.pupil_field).min()>=0.0 
        and np.real(coronagraph.pupil_field).max()<=1.0 
        and np.real(coronagraph.fpm_field).min()>=0.0
        and np.real(coronagraph.fpm_field).max()<=1.0):
        return 0.0
    else:  
        ln_prl = -np.inf
    ln_prr = ln_prior(design_parameters, input_pars)
    if not np.isfinite(ln_prr):
        return -np.inf 
    return lnprl + ln_prr + np.log(mrt) 
    

def run_emcee():
    with open('sleec_pars.yml', 'r') as f:
        input_pars = yml.full_load(f)
    print("Loaded input parameters from user configuration file")
    coronagraph = sli.import_pars()
    output0 = coronagraph(hci.Wavefront(coronagraph.pupil_field
                                        , coronagraph.wvl))
    print("Initialized coronagraph")
    sla.pupil_slepian_r(coefs=np.ones(input_pars['k'])
                        , corngrph=coronagraph, file_exists=False
                        , filename=input_pars['filename'], m=input_pars['m']
                        , n2d=input_pars['n2d'], k=input_pars['k'], L=170
                        , method='SE', orthgnlty=0, num_dense_pts=4096)
    print("Slepian file created")
    ini_positns = np.empty((input_pars['nwalkers'], input_pars['ndim']))
    for i in range(input_pars['k']):
        par_range = 1e-3*(input_pars['optmz_pupil_coefs'][2]
                          - input_pars['optmz_pupil_coefs'][1])
        ini_positns[:, i] = ( 
                             input_pars['optmz_pupil_coefs'][0] 
                             + par_range
                               *np.random.randn(input_pars['nwalkers'])                             )
    ini_positns[:, input_pars['k']] = ( 
        input_pars['optmz_fpm_dia'][0] 
        + 1e-3*(input_pars['optmz_fpm_dia'][2]-input_pars['optmz_fpm_dia'][1])
          *np.random.randn(input_pars['nwalkers'])
                                   ) 
    ini_positns[:, input_pars['k']+1] = ( 
        input_pars['optmz_lyot_dia'][0] 
        + 1e-3*(input_pars['optmz_lyot_dia'][2]-input_pars['optmz_lyot_dia'][1])
          *np.random.randn(input_pars['nwalkers'])
                                   ) 
    print("Initial parameter values (walker positions) set")
    sampler = emc.EnsembleSampler(input_pars['nwalkers'], input_pars['ndim']
                                    , ln_probability, args=(coronagraph
                                                          , input_pars))
    print("Ensemble Sampler instantiated.  Will start optimization run now.")
    sampler.run_mcmc(ini_positns, input_pars['nsteps'], progress=True)


if __name__ == '__main__':  
    run_emcee()






