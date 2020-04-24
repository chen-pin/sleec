"""<sleec_run>:  Module for running optimizations."""

import numpy as np
import emcee as emc
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
if (optmz_pupil_coefs[1] < coefs.all < optmz_pupil_coefs[2] 
    and optmz_fpm_dia[1] < fpm_dia < optmz_fpm_dia[2]
    and optmz_lyot_dia[1] < lyot_dia < optmz_lyot_dia[2]):
    return 0.0
return -np.inf
 

def ln_probability(design_parameters, coronagraph, input_pars):
"""
Return log probability for optimization using EMCEE.  

"""
    sla.pupil_slepian_r(coefs=design_parameters[0], corngraph=coronagraph
                        , file_exists=True, filename=input_pars['filename']
                        , m=None, n2d=None, k=None, L=None, method=None
                        , orthgnlty=None, num_dense_pts=None)
    sla.fpm_dot(fpm_dia=design_parameters[1], corngrph=coronagraph
                , num_airy=input_pars['num_airy'], field_stop=True)
    sla.lyot_circ_clear(lyot_dia=design_parameters[2])
    off_axis = (input_pars['IWA'] + input_pars['OWA'])/2.0
    mrt = slm.merit(coronagraph=coronagraph, inner_radius=input_pars['IWA']
                      , outer_radius=input_pars['OWA'], amplitude1=1.0
                      , area_penalty=input_pars['area_penalty']
                      , azimuthal_locatn=2*np.pi/7, amplitude2=1.0
                      , flux_ratio=input_pars['flux_ratio'])
    ln_prr = ln_prior(design_parameters, input_pars)
    if not np.isfinite(ln_prr):
        return -np.inf 
    return ln_prr + np.log(mrt) 
    

    
