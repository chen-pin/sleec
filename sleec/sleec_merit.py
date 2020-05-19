"""<sleec_merit>.  Module to compute merit/cost functions."""

import time
import numpy as np  
from scipy.special import jv
import hcipy as hci


def wvfrnt_annls_cnj(inner_radius, outer_radius, amplitude
                     , pupil_dia, grid, wvl):
    """
    Return a pupil wavefront that generates an annulus in the image.  
    
    Parameters
    ----------
    inner_radius : float
        Inner radius of the annulus in units of wvl/pupil_dia
    outer_radius : float
        Outer radius of the annulus in units of wvl/pupil_dia 
    amplitude : float
        Field amplitude factor
    pupil_dia : float
        Diameter of the pupil 
    grid : `hcipy.grid`  
        Pupil grid
    wvl : float
        Wavelength of light

    Returns
    -------
    return : hcipy.Wavefront
        Modulated wavefront illuminating the annuluar area in the image plane.

    """ 
    polar_grid = grid.as_('polar')
    argument = 2*np.pi*(outer_radius/pupil_dia)*polar_grid.r
    field1 = (outer_radius**2)*jv(1, argument)/argument 
    argument = 2*np.pi*(inner_radius/pupil_dia)*polar_grid.r
    field2 = (inner_radius**2)*jv(1, argument)/argument
    return hci.Wavefront(hci.Field(amplitude*(field1-field2), grid), wvl) 


def wvfrnt_off_axis(off_axis_locatn, azimuthal_locatn, amplitude, pupil_dia
                   , grid, wvl): 
    """
    Return a pupil wavefront to generate an off-axis PSF.  

    Parameters
    ----------
    off_axis_locatn : float
        Off-axis displacement in the image in units of wvl/pupil_dia
    azimuthal_locatn : float
        Azimuthal location in radians
    amplitude : float
        Field amplitude factor
    pupil_dia : float
        Pupil diameter
    grid : 'hcipy.grid'
        Pupil grid
    wvl : float
        Wavelength of light

    Returns
    -------
    return : hcipy.Wavefront 
        Modulated wavefront producing an off-axis PSF in the image plane 
    
    """ 
    polar_grid = grid.as_('polar')
    field = hci.Field(amplitude*np.exp(-1j*2*np.pi*(polar_grid.r/pupil_dia)
                             *off_axis_locatn
                             *np.cos(azimuthal_locatn-polar_grid.theta)), grid)
    return hci.Wavefront(field, wvl)


def throughput_annul(coronagraph, inner_radius, outer_radius, amplitude1
               , throughputs=False, insert_planet=True, **kwargs):
    """
    Return throughput in an annular area in the image.  

    Parameters:
    ----------
    coronagraph :  `sleec_optics.Coronagraph`
        The coronagraph obejct for optical propagation
    inner_radius : float
        Inner radius of the annulus in units of wvl/pupil_dia
    outer_radius : float
        Outer radius of the annulus in units of wvl/pupil_dia
    amplitude1 : float  
        Amplitude factor of annular illumination
    insert_planet : boolean 
        If `True` an offset PSF is generated in the image, and the keyword 
        arguments below are required
    throughputs : boolean
        If `True`, compute and return throughputs in each airy band.
    off_axis_locatn : float
        Off-axis displacement in the image in units of wvl/pupil_dia
    azimuthal_locatn : float
        Azimuthal location in radians
    amplitude2 : float
        Field amplitude factor

    Returns:
    -------
    whole_throughput : float
        Throughput averaged over the entire annulus  
    band_throughputs : ndarray 
        An array of throughput values.  Each averaged over an annular band 
        1 lambda/D in width, starting from the inner radius. If 
        (outer_radius - inner_radius) is not a whole number, the last band is
        a fractional band; only part of the band is illuminated but the 
        averaging is over the entire 1 lambda/D annulus. 
    planet_core_area : scalar  
        Number of pixels in the half-maximum encircled area of the planet's 
        final PSF 

    """
    # Synthesize sum pupil wavefront to produce annular illumination and an 
    # off-axis PSF in the focal plane
    field = ( 
                  (wvfrnt_annls_cnj(inner_radius, outer_radius, amplitude1
                                   , coronagraph.pupil_diamtr
                                   , coronagraph.pupil_field.grid
                                   , coronagraph.wvl)).electric_field 
            ) 
    if insert_planet:  
        field +=  (           
                     (wvfrnt_off_axis(kwargs['off_axis_locatn']
                                      , kwargs['azimuthal_locatn']
                                      , kwargs['amplitude2']
                                      , coronagraph.pupil_diamtr
                                      , coronagraph.pupil_field.grid
                                      , coronagraph.wvl)).electric_field
                  ) 
    wavefront = hci.Wavefront(field, coronagraph.wvl)
    output = coronagraph(wavefront)
    # Compute throughputs
    radii = output['on_fpm'].electric_field.grid.as_('polar').r  
    bounds = ( np.arange(inner_radius, outer_radius+1)
              *coronagraph.f_num*coronagraph.wvl )
    whole_denominator = np.sum(output['on_fpm'].intensity[np.where(
                                       (inner_radius*coronagraph.f_num
                                       *coronagraph.wvl <= radii)
                                       & (radii < outer_radius*coronagraph.f_num
                                                  *coronagraph.wvl) 
                                                       )])
    if (np.isfinite(whole_denominator) and whole_denominator != 0.0):
        pass
    else:
        print(
          "Invalid `whole_denominator` value, set it to `np.finfo(float).eps`"
             )
        whole_denominator = np.finfo(float).eps
    whole_numerator = np.sum(output['camera'].intensity[np.where(
                                       (inner_radius*coronagraph.f_num
                                       *coronagraph.wvl <= radii)
                                       & (radii < outer_radius*coronagraph.f_num
                                                  *coronagraph.wvl) 
                                                       )])
    whole_throughput = whole_numerator/whole_denominator
    if insert_planet:
        planet_intensity = output['camera'].intensity
        planet_core_area = np.where(planet_intensity>=planet_intensity.max()/2.0
                                    , 1, 0).sum()
    else:  
        planet_core_area = None 
    if throughputs:
        throughputs = [] 
        for i in range(len(bounds)-1):
            denominator = np.sum(output['on_fpm'].intensity[
                                                   np.where(
                                                       (bounds[i] <= radii) 
                                                        & (radii < bounds[i+1])
                                                           )
                                                           ]
                                )
            numerator = np.sum(
                               output['camera'].intensity[
                                                   np.where(
                                                       (bounds[i] <= radii) 
                                                        & (radii < bounds[i+1])
                                                           ) 
                                                         ]
                              )
            if (np.isfinite(denominator) and denominator != 0.0):
                pass
            else:
                print(
                  "Invalid `denominator` value, set it to `np.finfo(float).eps`"
                     )
                denominator = np.finfo(float).eps
            throughputs.append(numerator/denominator)
        return whole_throughput, np.array(throughputs), planet_core_area
    else:
        return whole_throughput, planet_core_area


def throughput_rand(coronagraph, inner_radius, outer_radius):
    """
    Return throughput for a planet in an annular area in the image.  The 
    planet's image position is randomly generated within the annular zone 
    between the inner_radius and outer_radius.   

    Parameters
    ----------
    coronagraph :  `sleec_optics.Coronagraph`
        The coronagraph obejct for optical propagation
    inner_radius : float
        Inner radius of the annulus in airy units (wvl/pupil_dia)
    outer_radius : float
        Outer radius of the annulus in airy units (wvl/pupil_dia)

    Returns
    -------
    throughput : float
        Throughput  
    planet_core_area : scalar  
        Number of pixels in the half-maximum encircled area of the planet's 
        final PSF 

    Notes
    -----
    - Throughput here is defined as the ratio of the intensity in the final 
      'camera' plane to the intensity incident on the first focal plane 
      (fpm).  Both intensities are integrated over the half-maximum core area 
      defined in the *camera plane*.  Therefore, decreasing the Lyot-stop 
      diameter incurrs a double penalty, from the smaller aperture area AND 
      from spreading out the planet image in the camera plane vs. that in 
      the fpm plane.  
    """
    r = inner_radius + (outer_radius-inner_radius)*np.random.random()
    theta = 2.0*np.pi*np.random.random()
    field =  (           
                 (wvfrnt_off_axis(off_axis_locatn=r
                                  , azimuthal_locatn=theta
                                  , amplitude=1.0
                                  , pupil_dia=coronagraph.pupil_diamtr
                                  , grid=coronagraph.pupil_field.grid
                                  , wvl=coronagraph.wvl)).electric_field
             ) 
    wavefront = hci.Wavefront(field, coronagraph.wvl)
    output = coronagraph(wavefront)
    # Compute throughputs
#    radii = output['on_fpm'].electric_field.grid.as_('polar').r  
#    bounds = ( np.arange(inner_radius, outer_radius+1)
#              *coronagraph.f_num*coronagraph.wvl )
#    throughputs = [] 
#    for i in range(len(bounds)-1):
#        denominator = np.sum(output['on_fpm'].intensity[
#                                               np.where((bounds[i] <= radii) 
#                                                       & (radii < bounds[i+1])) 
#                                              ])
#        numerator = np.sum(output['camera'].intensity[
#                                               np.where((bounds[i] <= radii) 
#                                                       & (radii < bounds[i+1])) 
#                                              ])
#        if (np.isfinite(denominator) and denominator != 0.0):
#            pass
#        else:
#            print(
#                  "Invalid `denominator` value, set it to `np.finfo(float).eps`"
#                 )
#            denominator = np.finfo(float).eps
#        throughputs.append(numerator/denominator)
#    whole_denominator = np.sum(output['on_fpm'].intensity[np.where(
#                                       (inner_radius*coronagraph.f_num
#                                       *coronagraph.wvl <= radii)
#                                       & (radii < outer_radius*coronagraph.f_num
#                                                  *coronagraph.wvl) 
#                                                       )])

    camera_intensity = output['camera'].intensity
    planet_core_ndx = np.where(camera_intensity>=camera_intensity.max()/2.0)
    planet_core_area = np.zeros_like(camera_intensity)
    planet_core_area[planet_core_ndx] = 1
    planet_core_area = planet_core_area.sum()
    numerator = np.sum(camera_intensity[planet_core_ndx])
    on_fpm_intensity = output['on_fpm'].intensity
    denominator = np.sum(on_fpm_intensity[planet_core_ndx])
#    if (np.isfinite(denominator) and denominator != 0.0):
#        pass
#    else:
#        print(
#          "Invalid `denominator` value, set it to `np.finfo(float).eps`"
#             )
#        denominator = np.finfo(float).eps
    throughput = numerator/denominator

#    import matplotlib.pyplot as plt
#    subplot_seq = 1
#    plt.figure(figsize=(15, 5))
#    for key in coronagraph.names:
#        print(key)
#        wf = output[key]
#        print(wf)
#        plt.subplot(2, len(coronagraph.names)//2, subplot_seq)
#        hci.imshow_field(wf.intensity)
#        plt.colorbar()
#        subplot_seq += 1
#    plt.show()
#    plt.close()

    return throughput, planet_core_area


def contrast(coronagraph, inner_radius, outer_radius, area_penalty
             , output_wavefronts, planet_core_area=None):
    """
    Return contrast in an annular area in the image.  

    Parameters:
    ----------
    coronagraph :  `sleec_optics.Coronagraph`
        The coronagraph obejct for optical propagation
    inner_radius : float
        Inner radius of the annulus in units of wvl/pupil_dia
    outer_radius : float
        Outer radius of the annulus in units of wvl/pupil_dia
    area_penalty : boolean
        If `True`, multiplies the raw contrast by the ratio of the number 
        pixels in the planet's half-intensity core to the number of pixels in 
        the stellar PSF core. 
    out_wavefronts : boolean
        If `True`, return propagated wavefronts.   
    planet_core_area : scalar  
        If `area_penalty` is `True`, this is th number of pixels in the 
        planet's half-intensity core (in the final image plane).  Otherwise, 
        it can be `None`.  

    Returns:
    -------
    whole_contrast : float
        Contrast averaged over the entire annulus  
    band_throughputs : ndarray 
        An array of contrast values.  Each averaged over an annular band 
        1 lambda/D in width, starting from the inner radius. If 
        (outer_radius - inner_radius) is not a whole number, the last band is
        a fractional band; only part of the band is illuminated but the 
        averaging is over the entire 1 lambda/D annulus. 
    penalty : float
        The penalty factor associated with `area_penalty`.  

    """
    output = coronagraph(hci.Wavefront(coronagraph.pupil_field
                                       , coronagraph.wvl))
    psf = output['on_fpm'].intensity
    if (np.isfinite(psf.max()) and psf.max()!=0.0):
        psf0 = psf.max()
    else:
        print("Invalid `psf0` value, set it to `np.finfo(float).eps`")
        psf0 = np.finfo(float).eps
    radii = output['on_fpm'].electric_field.grid.as_('polar').r  
    bounds = ( np.arange(inner_radius, outer_radius+1)
              *coronagraph.f_num*coronagraph.wvl )
    contrasts = [] 
    for i in range(len(bounds)-1):
        mask = np.where((bounds[i] <= radii) & (radii < bounds[i+1]), 1, 0)
        avg_background = ( np.sum(output['camera'].intensity[np.where(mask)])
                       /mask.sum() )
        contrasts.append(avg_background/psf0)
    mask = np.where((inner_radius*coronagraph.f_num*coronagraph.wvl <= radii)
                    & (radii < outer_radius*coronagraph.f_num*coronagraph.wvl)
                    , 1, 0)
    whole_background = ( np.sum(output['camera'].intensity[np.where(mask)])
                       /mask.sum() )
    whole_contrast = whole_background/psf0 
    if area_penalty:
        psf_core_area = np.where(psf>=psf0/2.0, 1, 0).sum()
        if (np.isfinite(psf_core_area) and psf_core_area != 0):
            penalty = planet_core_area/psf_core_area
        else:
            print("Invalid `psf0' value, set `penalty` to 1E15")
            penalty = 1e15
        whole_contrast *= penalty
    if output_wavefronts:
        return whole_contrast, np.array(contrasts), penalty, output
    else:
        return whole_contrast, np.array(contrasts), penalty


def merit(coronagraph, inner_radius, outer_radius, amplitude1
         , area_penalty, contrast_weights, off_axis_locatn, azimuthal_locatn
         , amplitude2, flux_ratio, output_wavefronts):
    """
    Return a merit value, throughputs, and contrasts

    Parameters:
    ----------
    coronagraph :  `sleec_optics.Coronagraph`
        The coronagraph obejct for optical propagation
    inner_radius : float
        Inner radius of the annulus in units of wvl/pupil_dia
    outer_radius : float
        Outer radius of the annulus in units of wvl/pupil_dia
    amplitude1 : float  
        Amplitude factor of annular illumination
    area_penalty : boolean
        If `True`, multiplies the raw contrast by the ratio of the number 
        pixels in the planet's half-intensity core to the number of pixels in 
        the stellar PSF core.  
    contrast_weights : `None` or array_like
        Weights to be applied to each airy band of contrast.  If `None`, 
        use `whole_contrast`.  
    off_axis_locatn : float
        Off-axis displacement in the image in units of wvl/pupil_dia
    azimuthal_locatn : float
        Azimuthal location in radians
    amplitude2 : float
        Field amplitude factor
    flux_ratio : float
        Planet-to-star flux ratio
    out_wavefronts : boolean
        If `True`, return propagated wavefronts.   

    Returns
    -------
    merit : float
        Merit value
    throughput_output : tuple
        Outputs from the `throughput()` function
    contrast_output : tuple
        Outputs from the `contrast()` function 

    """

    throughput_output = throughput(coronagraph, inner_radius, outer_radius)
    contrast_output = contrast(coronagraph, inner_radius, outer_radius
                               , area_penalty
                               , output_wavefronts=output_wavefronts
                               , planet_core_area=throughput_output[1])
    eta = throughput_output[0]
    if contrast_weights == None:
        C = contrast_output[0]
    else:
        C = (contrast_output[1]*contrast_weights).mean()
    epsilon = flux_ratio
    if (np.isfinite(eta) and np.isfinite(C) and eta!=0.0 and C!=0.0):
        merit = eta/np.sqrt(eta+C/epsilon)
    else:
        print("Invalid values for computing `merit`, set it to machine eps")
        merit = np.finfo(float).eps
    if output_wavefronts:
        return merit, throughput_output, contrast_output, contrast_output[3]
    print("eta = {}\n C= {}\n epsilon = {}\n merit= {}".format(eta, C
                                                                 , epsilon
                                                                 , merit))
    return merit, throughput_output, contrast_output


def _test_wvfrnt_annls_cnj():
    import matplotlib.pyplot as plt
    import yaml as yml
    from matplotlib.colors import ListedColormap
    from turbo_colormap import turbo_colormap_data
    import sleec_input as sli

    turbo = ListedColormap(turbo_colormap_data)

    coronagraph, pars = sli.import_pars()
    print('Coronagraph initialized')
    output = coronagraph(wvfrnt_annls_cnj(pars['IWA'], pars['OWA'], 1.0
                         , coronagraph.pupil_diamtr
                         , coronagraph.pupil_field.grid
                         , coronagraph.wvl))  
    print('Propagation complete')
    subplot_seq = 1
    plt.figure(figsize=(15, 5))
    for key in coronagraph.names:
        print(key)
        wf = output[key]
        print(wf)
        plt.subplot(2, len(coronagraph.names)//2, subplot_seq)
        hci.imshow_field(wf.intensity, cmap=turbo)
        plt.colorbar()
        subplot_seq += 1
    plt.show()
    plt.close()


def _test_annls_offaxis():
    import matplotlib.pyplot as plt
    import yaml as yml
    from matplotlib.colors import ListedColormap
    from turbo_colormap import turbo_colormap_data
    import sleec_input as sli

    turbo = ListedColormap(turbo_colormap_data)

    coronagraph, pars = sli.import_pars()
    print('Coronagraph initialized')
    print('Coronagraph initialized')
    field = (
                 (wvfrnt_annls_cnj(pars['IWA'], pars['OWA'], 1.0
                         , coronagraph.pupil_diamtr
                         , coronagraph.pupil_field.grid
                         , coronagraph.wvl)).electric_field  
                + (wvfrnt_off_axis(3.0, np.pi/4.0, 1.0
                         , coronagraph.pupil_diamtr
                         , coronagraph.pupil_field.grid
                         , coronagraph.wvl)).electric_field
                )  
    wavefront = hci.Wavefront(field, coronagraph.wvl)
    output = coronagraph(wavefront)
    print('Propagation complete')
    subplot_seq = 1
    plt.figure(figsize=(15, 5))
    for key in coronagraph.names:
        print(key)
        wf = output[key]
        print(wf)
        plt.subplot(2, len(coronagraph.names)//2, subplot_seq)
        hci.imshow_field(wf.intensity, cmap=turbo)
        plt.colorbar()
        subplot_seq += 1
    plt.show()
    plt.close()


def _test_throughput():
    import sleec_input as sli
    coronagraph, pars = sli.import_pars()
    print('Coronagraph initialized')
    output = throughput(coronagraph, pars['IWA'], pars['OWA'])
    print("Average throughput over entire annulus:  {}".format(output[0]))
    print("Number of pixels in planet core:  {}".format(output[1]))


def _test_contrast():
    import yaml as yml
    import sleec_input as sli
    inner_radius = float(input('Enter inner radius in units of lambda/D:  '))
    outer_radius = float(input('Enter outer radius in units of lambda/D:  '))
    coronagraph, pars = sli.import_pars()
    print('Coronagraph initialized')
    output = contrast(coronagraph, inner_radius, outer_radius, True, 54)
    print("Average contrast over entire annulus:  {}".format(output[0]))
    print("Core-footprint penalty:  {}".format(output[2]))
    print("Contrast in 1 lambda/D bands")
    for i in range(len(output[1])):
        print('Band {}:  {}'.format(i+1, output[1][i]))


def _test_merit():
    import yaml as yml
    import sleec_input as sli
    inner_radius = float(input('Enter inner radius in units of lambda/D:  '))
    outer_radius = float(input('Enter outer radius in units of lambda/D:  '))
    coronagraph, pars = sli.import_pars()
    print('Coronagraph initialized')
    off_axis_locatn = (inner_radius+outer_radius)/2.0
    m = merit(coronagraph, inner_radius, outer_radius, amplitude1=1.0
         , area_penalty=True, off_axis_locatn=off_axis_locatn
         , azimuthal_locatn=2*np.pi/7
         , amplitude2=1, flux_ratio=1e-10)
    print("Merit value = {}".format(m[0]))
    print("Throughput = {}".format(m[1][0]))
    print("Contrast = {}".format(m[2][0]))


if __name__ == '__main__': 
#    _test_wvfrnt_annls_cnj()
    _test_throughput()
#    _test_annls_offaxis()
#    _test_contrast()
#    _test_merit()
