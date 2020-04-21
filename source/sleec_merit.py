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


def throughput(coronagraph, inner_radius, outer_radius, pupil_dia):
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
    pupil_dia : float
        Diameter of the pupil 

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

    Notes
    -----
    - Need to implement propagation with clear pupil to compute denominators.  

    """
    field = (
                 (wvfrnt_annls_cnj(pars['IWA'], pars['OWA'], 1.0
                         , pars['pupil_dia'], coronagraph.pupil_field.grid
                         , coronagraph.wvl)).electric_field  
                + (wvfrnt_off_axis((pars['IWA']+pars['OWA'])/2.0, np.pi/4.0, 1.0
                         , pars['pupil_dia'], coronagraph.pupil_field.grid
                         , coronagraph.wvl)).electric_field
                )  
    wavefront = hci.Wavefront(field, coronagraph.wvl)
    output = coronagraph(wvfrnt_annls_cnj(inner_radius, outer_radius
                                          , pupil_dia
                                          , coronagraph.pupil_field.grid
                                          , coronagraph.wvl) )
    radii = output['on_fpm'].electric_field.grid.as_('polar').r  
    bounds = ( np.arange(inner_radius, outer_radius+1)
              *coronagraph.f_num*coronagraph.wvl )
    throughputs = [] 
    for i in range(len(bounds)-1):
        denominator = np.sum(output['on_fpm'].intensity[
                                               np.where((bounds[i] <= radii) 
                                                       & (radii < bounds[i+1])) 
                                              ])
        numerator = np.sum(output['camera'].intensity[
                                               np.where((bounds[i] <= radii) 
                                                       & (radii < bounds[i+1])) 
                                              ])
        throughputs.append(numerator/denominator)
    whole_denominator = np.sum(output['on_fpm'].intensity[np.where(
                                       (inner_radius*coronagraph.f_num
                                       *coronagraph.wvl <= radii)
                                       & (radii < outer_radius*coronagraph.f_num
                                                  *coronagraph.wvl) 
                                                       )])
    whole_numerator = np.sum(output['camera'].intensity[np.where(
                                       (inner_radius*coronagraph.f_num
                                       *coronagraph.wvl <= radii)
                                       & (radii < outer_radius*coronagraph.f_num
                                                  *coronagraph.wvl) 
                                                       )])
    whole_throughput = whole_numerator/whole_denominator
    return whole_throughput, np.array(throughputs)


def _test_wvfrnt_annls_cnj():
    import matplotlib.pyplot as plt
    import yaml as yml
    from matplotlib.colors import ListedColormap
    from turbo_colormap import turbo_colormap_data
    import sleec_input as sli

    turbo = ListedColormap(turbo_colormap_data)

    with open('sleec_pars.yml', 'r') as f:
        pars = yml.full_load(f)
    print('Loaded parameters from .yml file')
    coronagraph = sli.import_pars()
    print('Coronagraph initialized')
    output = coronagraph(wvfrnt_annls_cnj(pars['IWA'], pars['OWA'], 1.0
                         , pars['pupil_dia'], coronagraph.pupil_field.grid
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

    with open('sleec_pars.yml', 'r') as f:
        pars = yml.full_load(f)
    print('Loaded parameters from .yml file')
    coronagraph = sli.import_pars()
    print('Coronagraph initialized')
    field = (
                 (wvfrnt_annls_cnj(pars['IWA'], pars['OWA'], 1.0
                         , pars['pupil_dia'], coronagraph.pupil_field.grid
                         , coronagraph.wvl)).electric_field  
                + (wvfrnt_off_axis(3.0, np.pi/4.0, 1.0
                         , pars['pupil_dia'], coronagraph.pupil_field.grid
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
    import yaml as yml
    import sleec_input as sli
    inner_radius = float(input('Enter inner radius in units of lambda/D:  '))
    outer_radius = float(input('Enter outer radius in units of lambda/D:  '))
    with open('sleec_pars.yml', 'r') as f:
        pars = yml.full_load(f)
    print('Loaded parameters from .yml file')
    coronagraph = sli.import_pars()
    print('Coronagraph initialized')
    output = throughput(coronagraph, inner_radius, outer_radius
                        , pars['pupil_dia']
                       ) 
    print("Average throughput over entire annulus:  {}".format(output[0]))
    print("Throughput in 1 lambda/D bands")
    for i in range(len(output[1])):
        print('Band {}:  {}'.format(i+1, output[1][i]))


if __name__ == '__main__': 
#    _test_wvfrnt_annls_cnj()
#    _test_throughput()
    _test_annls_offaxis()
