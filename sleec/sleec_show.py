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
    contrast_map = (wavefronts['camera'].intensity
                    /wavefronts['on_fpm'].intensity.max())
    print("merit value:  {}".format(merit_val))
    print("required exposure time:  {}".format(tau))
    print("overall throughput:  {}".format(whole_throughput))
    print("throughput each Airy band:  {}".format(throughputs))
    print("overall contrast:  {}".format(whole_contrast))
    print("contrast in each Airy band:  {}".format(band_contrasts))
    print("encircled PSF energy {}".format(encircle_psf))
    grid_scale = input_pars['f_num']*input_pars['wvl']
    hci.imshow_field(coronagraph.pupil_field, cmap=turbo
                     , grid_units=input_pars['pupil_dia'])
    plt.title('PUPIL FIELD')
    plt.colorbar()
    plt.show()
    hci.imshow_field(coronagraph.fpm_field, cmap=turbo, grid_units=grid_scale)
    plt.title('FPM FIELD')
    plt.colorbar()
    plt.show()
    hci.imshow_field(coronagraph.lyot_field, cmap=turbo
                     , grid_units=input_pars['pupil_dia'])
    plt.title('LYOT FIELD')
    plt.colorbar()
    plt.show()
    hci.imshow_field(np.log10(contrast_map), cmap=turbo
                     ,vmax=0.0, vmin=-15.0, grid_units=grid_scale)
    plt.title('LOG CONTRAST')
    plt.colorbar()
    plt.show()


if __name__ == "__main__":
#    show_run_stats()    
#    reconstruct('20200508t1759', 35937)
#    reconstruct('20200516t0007', 10712)
#    reconstruct('20200516t0007', 10732)
#    reconstruct('20200516t0007', 27391)
#    reconstruct('20200518t0951', 7211)
    reconstruct('20200706t2300', 12922)
