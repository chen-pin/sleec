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
    burnin = int(2 * np.max(tau))
    thin = int(0.5 * np.min(tau))
#    samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
    samples = reader.get_chain(discard=burnin, flat=True)
#    log_prob_samples = reader.get_log_prob(discard=burnin, flat=True, thin=thin)
    log_prob_samples = reader.get_log_prob(discard=burnin, flat=True)
    log_prob_max = log_prob_samples.max()
    log_prob_max_ndx = np.argmax(log_prob_samples)
    print("auto-correlation times:  {}".format(tau))
    print("burn-in: {0}".format(burnin))
    print("thin: {0}".format(thin))
    print("flat chain shape: {0}".format(samples.shape))
    print("flat log prob shape: {0}".format(log_prob_samples.shape))
    print("maximum log prob value:  {}".format(log_prob_max))
    print("index of max value:  {}\n".format(log_prob_max_ndx))
    ndim = samples.shape[-1]
    blob = reader.get_blobs(discard=burnin, flat=True)
    print("blob shape:  {}: ".format(blob.shape))
    sortarg = np.flip(np.argsort(log_prob_samples))
    print("TOP FIVE SAMPLES")
    for j in range(5):
        print("sample #  {}".format(sortarg[j]))
        print("ln probability: {} ".format(log_prob_samples[sortarg[j]]))
        print("throughput: {} ".format(blob['whole_throughput'][sortarg[j]]))
        print("contrast: {} ".format(blob['whole_contrast'][sortarg[j]]))
        print("planet core area: {} ".format(blob['planet_core_area']
                                               [sortarg[j]]))
        print("area penalty: {} ".format(blob['area_penalty'][sortarg[j]]))
        for i in range(ndim): 
            print("parameter {} value:  {}"
                  .format(i, samples[sortarg[j], i]))
    all_samples = np.concatenate(
        (samples, log_prob_samples[:, None]), axis=1
                                )
    labels = list(map(r"$\theta_{{{0}}}$".format, range(1, ndim + 1)))
    labels += ["log prob"]
    corner.corner(all_samples, labels=labels)
    plt.show()


def reconstruct(filename, sample_number):  
    path = os.path.join('..', '..', 'sleec_out', filename)
    reader = emc.backends.HDFBackend(path)
    pars = reader.get_chain(flat=True)
    with open('sleec_pars.yml', 'r') as f:
        input_pars = yml.full_load(f)
    coronagraph = sli.import_pars()
    (ln_prob, whole_throughput, planet_core_area, whole_contrast, area_penalty
     , coronagraph, wavefronts, throughputs, contrasts) \
     = slr.ln_probability(design_parameters=pars, coronagraph=coronagraph
                          , input_pars=input_pars, reconstruct=True) 
    contrast_map = (wavefronts['camera'].intensity
                    /wavefronts['on_fpm'].intensity.max())
    print("ln probability:  {}".format(ln_prob))
    print("overall throughput:  {}".format(whole_throughput))
    print("# of pixels in planet core:  {}".format(planet_core_area))
    print("(planet cors)/(PSF core):  {}".format(whole_contrast))
    print("overall contrast:  {}".format(whole_contrast))
    print("throughput in each airy band {}".format(throughputs))
    print("contrast in each airy band {}".format(contrasts))
    hci.imshow_field(coronagraph.pupil_field, cmap=turbo)
    plt.title('PUPIL FIELD')
    plt.colorbar()
    plt.show()
    hci.imshow_field(coronagraph.fpm_field, cmap=turbo)
    plt.title('FPM FIELD')
    plt.colorbar()
    plt.show()
    hci.imshow_field(coronagraph.lyot_field, cmap=turbo)
    plt.title('LYOT FIELD')
    plt.colorbar()
    plt.show()
    hci.imshow_field(contrast_map, cmap=turbo)
    plt.title('CONTRAST')
    plt.colorbar()
    plt.show()
    

if __name__ == "__main__":
    show_run_stats()    
