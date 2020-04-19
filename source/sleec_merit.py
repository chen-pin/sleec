import numpy as np  
from scipy.special import jv
import hcipy as hci


def wvfrnt_annls_cnj(inner_radius, outer_radius, pupil_dia, grid, wvl):
    """
    Return a pupil wavefront that generates an annulus in the image.  
    
    Parameters
    ----------
    inner_radius : float
        Inner radius of the annulus in units of wvl/pupil_dia
    outer_radius : float
        Outer radius of the annulus in units of wvl/pupil_dia 
    pupil_dia : float
        Diameter of the pupil 
    grid : `hcipy.grid`  
        Pupil grid
    wvl : float
        Wavelength of light

    """ 
    polar_grid = grid.as_('polar')
    argument = 2*np.pi*(outer_radius/pupil_dia)*polar_grid.r
    field1 = (outer_radius**2)*jv(1, argument)/argument 
    argument = 2*np.pi*(inner_radius/pupil_dia)*polar_grid.r
    field2 = (inner_radius**2)*jv(1, argument)/argument
    return hci.Wavefront(hci.Field((field1-field2), grid), wvl) 


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
    output = coronagraph(wvfrnt_annls_cnj(pars['IWA'], pars['OWA']
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


if __name__ == '__main__': 
    _test_wvfrnt_annls_cnj()
