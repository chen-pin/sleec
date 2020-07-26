import os
import yaml as yml
import numpy as np
import matplotlib.pyplot as plt
import corner
import emcee as emc
import hcipy as hci
import sleec_run as slr
import sleec_input as sli
from matplotlib.colors import ListedColormap
from turbo_colormap import turbo_colormap_data
turbo = ListedColormap(turbo_colormap_data)


def show_run_stats():
    filename = input('Enter EMCEE backend filename:  ')
    path = os.path.join('..', '..', 'sleec_out', filename)
    reader = emc.backends.HDFBackend(path)
    tau = reader.get_autocorr_time(tol=4)
    samples = reader.get_chain(flat=True)
    log_prob_samples = reader.get_log_prob(flat=True)
    log_prob_max = log_prob_samples.max()
    log_prob_max_ndx = np.argmax(log_prob_samples)
    print("auto-correlation times:  {}".format(tau))
    print("flat chain shape: {0}".format(samples.shape))
    print("flat log prob shape: {0}".format(log_prob_samples.shape))
    print("maximum log prob value:  {}".format(log_prob_max))
    print("index of max value:  {}\n".format(log_prob_max_ndx))
    ndim = samples.shape[-1]
    blob = reader.get_blobs(flat=True)
    print("blob shape:  {}: ".format(blob.shape))
    sortarg = np.flip(np.argsort(log_prob_samples))
    ndx = [sortarg[0]]
    count = 0 
    prev_log_prob = log_prob_samples[sortarg[0]]
    for n in sortarg[1:]:
        log_prob = log_prob_samples[n]
        if (log_prob != prev_log_prob and count<=3):
            ndx.append(n)
            prev_log_prob = log_prob
            count += 1
        elif count>3:
            break
        else:
            continue
    print("\nTOP FIVE SAMPLES")
    for j in ndx:
        print("\nsample #  {}".format(j))
        print("ln probability: {} ".format(log_prob_samples[j]))
        print("throughput: {} ".format(blob['whole_throughput'][j]))
        print("contrast: {} ".format(blob['whole_contrast'][j]))
        print("encircled PSF energy: {} ".format(blob['encircle_psf']
                                               [j]))
        for i in range(ndim): 
            print("parameter {} value:  {}"
                  .format(i+1, samples[j, i]))
    all_samples = np.concatenate(
        (samples, log_prob_samples[:, None]), axis=1
                                )
    labels = list(map(r"$\theta_{{{0}}}$".format, range(1, ndim + 1)))
    labels += ["log prob"]
    corner.corner(all_samples, labels=labels)
    plt.show()
    save = input("Save figure? (y/[n]):  ")
    if save == 'y':
        filename = input("Enter filename:  ")
        path = os.path.join('..', '..', 'sleec_out', filename)
        plt.savefig(path)
    log_prob_samples = reader.get_log_prob(flat=False)
    iterations = log_prob_samples.shape[0]
    walkers = log_prob_samples.shape[1]
    x = [i for i in range(iterations)]
    for walker in range(walkers): 
        plt.plot(x, log_prob_samples[:, walker])
    plt.show()
    save = input("Save figure? (y/[n]):  ")
    if save == 'y':
        filename = input("Enter filename:  ")
        path = os.path.join('..', '..', 'sleec_out', filename)
        plt.savefig(path)


def quad_plot(field, suptitle, grid_scale, grid_label, scale1='linear'
              , im_vmin=None, im_vmax=None):
    intensity = (np.abs(field))**2
    if np.iscomplexobj(field):
        phase = np.angle(field)
    else:
        phase = np.zeros_like(field)
    plt.figure(figsize=(15, 12))
    plt.suptitle(suptitle)
    plt.subplot(221)
    plt.title('intensity')
    hci.imshow_psf(intensity, scale=scale1, normalization='none'
                   , cmap=turbo, vmin=im_vmin, vmax=im_vmax
                   , spatial_resolution=grid_scale)
    plt.subplot(223)
    plt.title('phase [pi]')
    hci.imshow_psf(phase/np.pi, scale='linear', normalization='none'
                   , cmap=turbo, vmax=1.0, vmin=-1.0
                   , spatial_resolution=grid_scale)
    plt.subplot(222)
    side_dim = int(np.sqrt(len(field)))
    ctr_ndx = side_dim//2
    plt.title('Horizontal Center Cross-Section:  Intensity')
    abscissa = ( field.grid.x.reshape(side_dim, side_dim)
                                                     [ctr_ndx, :] )/grid_scale
    plt.plot( abscissa 
              , intensity.reshape(side_dim, side_dim)[ctr_ndx,:] )
    plt.xlabel(grid_label)
    plt.ylabel('amplitude')
    plt.subplot(224)
    plt.title('Horizontal Center Cross-Section:  Phase')
    plt.plot( abscissa 
              , phase.reshape(side_dim, side_dim)[ctr_ndx,:]/np.pi )
    plt.xlabel(grid_label)
    plt.ylabel('phase [pi]')
    plt.show()
    
    

def reconstruct(fn_prefix, sample_number, airy_scale=True):  
    path = os.path.join('..', '..', 'sleec_out')
    file_path_hdf = os.path.join(path, fn_prefix+'.hdf')
    file_path_slp = os.path.join(path, fn_prefix+'.pkl')
    file_path_yml = os.path.join(path, fn_prefix+'.yml')
    reader = emc.backends.HDFBackend(file_path_hdf)
    pars = reader.get_chain(flat=True)[sample_number]
    coronagraph, input_pars = sli.import_pars(file_path_yml)
    (merit_val, tau, whole_throughput, throughputs, whole_contrast
     , band_contrasts, encircle_psf, wavefronts, Cr_p, Cr_b_star, Cr_b_zodi
     , Cr_b_exozodi, Cr_b) = slr.ln_probability(design_parameters=pars
                                 , coronagraph=coronagraph
                                 , non_negative=input_pars['non_negative']
                                 , input_pars=input_pars, reconstruct=True) 
    print("merit value:  {}".format(merit_val))
    print("required exposure time:  {}".format(tau))
    print("overall throughput:  {}".format(whole_throughput))
    print("throughput each Airy band:  {}".format(throughputs))
    print("overall contrast:  {}".format(whole_contrast))
    print("contrast in each Airy band:  {}".format(band_contrasts))
    print("encircled PSF energy {}".format(encircle_psf))
    # Check Contrast
    iwa = input_pars['IWA']
    owa = input_pars['OWA']
    resel = coronagraph.f_num*coronagraph.wvl
#    r = coronagraph.fpm_field.grid.as_('polar').r
#    ndx = np.where((iwa*resel <= r)*( r  < owa*resel))
#    print("\nndx sum:  {}".format(np.sum(np.ones_like(r)[ndx])))
#
#    field = ( np.abs(wavefronts['camera'].electric_field)
#              /np.abs(wavefronts['on_fpm'].electric_field.max()) )
#    denom = np.abs(wavefronts['on_fpm'].electric_field)
#    truth = np.where((denom>=0.0)*(np.isreal(denom)), 1, 0).reshape(
#                int(np.sqrt(451584)), int(np.sqrt(451584)))
#    print("truth sum:  {}".format(truth.sum()))
#    print("denominator:  {}:".format(
#                        np.abs(wavefronts['on_fpm'].electric_field.max())
#                                      ))
#    print("denominator^2:  {}:".format(
#                        np.abs(wavefronts['on_fpm'].electric_field.max())**2
#                                      ))
#    contrast_1 = (np.abs(field))**2
#    contrast_zone = contrast_1[ndx].mean()
#    print("contrast_1:  {}".format(contrast_1))
#    print("Mean contrast_1 in target zone:  {:E}".format(contrast_zone))
#
#    field = np.abs(wavefronts['camera'].electric_field)**2
#    norm = (np.abs(wavefronts['on_fpm'].electric_field)**2).max() 
#    print("\nnorm:  {}".format(norm))
#    contrast_2 = (field/norm)
#    contrast_zone = contrast_2[ndx].mean()
#    print("contrast_2:  {}".format(contrast_2))
#    print("Mean contrast_2 in target zone:  {:E}".format(contrast_zone))
#
    psf = wavefronts['on_fpm'].intensity
    psf0 = psf.max()
#    print("\npsf0:  {}".format(psf0))
    radii = wavefronts['on_fpm'].electric_field.grid.as_('polar').r
    mask = np.where((iwa*resel <= radii) & (radii < owa*resel), 1, 0)
#    print("mask sum:  {}".format(mask.sum()))
#    avg_background = ( np.sum(wavefronts['camera'].intensity[np.where(mask)])
#                       /mask.sum() )
    avg_background = wavefronts['camera'].intensity[np.where(mask)].mean()
    contrast_r = avg_background/psf0
    print("recomputed contrast:  {:E}".format(contrast_r))
    # Show PUPIL
    field = coronagraph.pupil_field
    grid_scale = coronagraph.pupil_diameter
    quad_plot(field, 'PUPIL', grid_scale, 'x [pupil diameter]') 
    # Show ON_FPM
    field = wavefronts['on_fpm'].electric_field
    grid_scale = coronagraph.f_num*coronagraph.wvl
    quad_plot(field, 'INCIDENT FIELD ON THE FPM', grid_scale, 'lambda/D')
    # Show FPM
    field = coronagraph.fpm_field
    grid_scale = coronagraph.f_num*coronagraph.wvl
    quad_plot(field, 'FPM', grid_scale, '[lambda/D]') 
    # Show ON_LYOT
    field = wavefronts['on_lyot'].electric_field
    grid_scale = coronagraph.pupil_diameter
    quad_plot(field, 'FIELD INCIDENT ON THE LYOT PUPIL', grid_scale
              , 'x [pupil diameter]') 
    # Show LYOT APERTURE
    field = coronagraph.lyot_field
    grid_scale = coronagraph.pupil_diameter
    quad_plot(field, 'LYOT STOP', grid_scale, 'x [pupil diameter]') 
    # Show CONTRAST	
    field = ( np.abs(wavefronts['camera'].electric_field)
              /np.abs(wavefronts['on_fpm'].electric_field).max() )
    scale_max = 10**(np.ceil(np.log10(field.max()**2)))
    scale_min = 1e-13
    grid_scale = coronagraph.f_num*coronagraph.wvl
    quad_plot(field, 'CONTRAST', grid_scale, '[lambda/D]', 'log'
              , im_vmin=scale_min, im_vmax=scale_max) 


if __name__ == "__main__":
    show_run_stats()    
#    reconstruct('20200508t1759', 35937)
#    reconstruct('20200516t0007', 10712)
#    reconstruct('20200516t0007', 10732)
#    reconstruct('20200516t0007', 27391)
#    reconstruct('20200518t0951', 7211)
#    reconstruct('20200706t2300', 12942)
