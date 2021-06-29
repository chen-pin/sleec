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


class RunData(emc.backends.HDFBackend):
    """
    Retrieve and plot EMCEE run results.  

    Parameters
    ----------
    run_id : str
        Run identifier (e.g. 20210624t1612)  
    read_only : bool, optional
        If `True`, the backend will throw a `RuntimeError` if the file is 
        opened with write access
    directory : str
        Base data-directory path

    Attributes
    ----------
    run_id : str
        Run identifier (e.g. 20210624t1612)  
    read_only : bool, optional
        If `True`, the backend will throw a `RuntimeError` if the file is 
        opened with write access
    directory : str
        Base data-directory path
    path : str
        Full path to the data record
    samples : ndarray
        The MCMC samples, array[..., num_walkers, num_dims] of sampled
        parameter values
    log_prob_samples : ndarray
        Array[..., num_walkers] of ln probabilities evaluated at the MCMC 
        samples
    blob : ndarray
        Array[..., num_walkers] of blobs (appended saved info.) for each sample
    num_walkers : int
        Number of walkers in the EMCEE run
    num_dims : int
        Number of dimensions (parameters) in the EMCEE run
    num_iters : float
        Number of iterations executed in the EMCEE run

    """
    def __init__(self
                 , run_id, read_only=True
                 , directory=os.path.join('..', '..', 'sleec_out')):
        self.run_id = run_id
        self.read_only = read_only
        self.directory = directory
        self.path = os.path.join(directory, self.run_id+'.hdf')
        super().__init__(filename=self.path, read_only=read_only)
        self.samples = self.get_chain(flat=True)
        self.log_prob_samples = self.get_log_prob(flat=True)
        self.blob = self.get_blobs(flat=True)
        self.num_walkers = self.shape[0]
        self.num_dims = self.shape[1]
        self.num_iters = self.samples.shape[0]/self.num_walkers

    def print_summary(self):
        """ 
        Print summary of EMCEE run characteristics
        """
        log_prob_max = self.log_prob_samples.max()
        ndx = np.argmax(self.log_prob_samples)
        print("number of walkers: {0}".format(self.num_walkers))
        print("number of dimensions: {0}".format(self.num_dims))
        print("number of iterations: {0}".format(self.num_iters))
        print("flat chain shape: {0}".format(self.samples.shape))
        print("flat log prob shape: {0}".format(self.log_prob_samples.shape))
        print("maximum log prob value:  {}".format(log_prob_max))
        print("index of max value:  {}".format(ndx))
        print("throughput: {} ".format(self.blob['whole_throughput'][ndx]))
        print("contrast: {} ".format(self.blob['whole_contrast'][ndx]))
        print("encircled PSF energy: {} "
              .format(self.blob['encircle_psf'][ndx]))
        print("auto-correlation times:  \n{}"
              .format(self.get_autocorr_time(tol=1)))

    def get_top_samples(self, n=5, print_summary=True, print_pars=False):
        """
        Return `self.samples` indices of the `n` samples with the highest 
        ln probability values.  

        """
        sortarg = np.flip(np.argsort(self.log_prob_samples))
        ndx = [sortarg[0]]
        count = 0 
        prev_log_prob = self.log_prob_samples[sortarg[0]]
        if print_summary:
            for i in sortarg[1:]:
                log_prob = self.log_prob_samples[i]
                if (log_prob != prev_log_prob and count<=n-2):
                    ndx.append(i)
                    prev_log_prob = log_prob
                    count += 1
                elif count>3:
                    break
                else:
                    continue
            print("\nSUMMARY OF THE TOP "+str(n)+" SAMPLES")
            for j in ndx:
                print("\nsample #  {}".format(j))
                print("ln probability: {} ".format(self.log_prob_samples[j]))
                print("throughput: {} "
                      .format(self.blob['whole_throughput'][j]))
                print("contrast: {} ".format(self.blob['whole_contrast'][j]))
                print("encircled PSF energy: {} "
                      .format(self.blob['encircle_psf']
                                                       [j]))
        if print_pars:
            for i in range(ndim): 
                print("parameter {} value:  {}"
                      .format(i+1, self.samples[j, i]))
        return ndx[:n]

    def corner_plot(self, save=False):
        all_samples = np.concatenate(
            (self.samples, self.log_prob_samples[:, None]), axis=1)
        labels = list(map(r"$\theta_{{{0}}}$".format
                      , range(1, self.num_dims + 1)))
        labels += ["log prob"]
        corner.corner(all_samples, labels=labels)
        if save:
            filename = self.run_id+'_corner.png'
            path = os.path.join(self.directory, filename)
            plt.savefig(path)
        else:
            plt.show()

    def track_plot(self, save=False):
        non_flat_log_probs = self.log_prob_samples.reshape(-1, self.shape[0])
        x = [i for i in range(int(self.num_iters))]
        for walker in range(self.num_walkers): 
            plt.plot(x, non_flat_log_probs[:, walker])
            plt.xlabel('iteration')
            plt.ylabel('ln probability')
        if save:
            filename = self.run_id+'_tracks.png'
            path = os.path.join(self.directory, filename)
            plt.savefig(path)
        else:
            plt.show()

    def reconstruct(self, sample_number, airy_scale=True):  
        file_path_slp = os.path.join(self.directory, self.run_id+'.pkl')
        file_path_yml = os.path.join(self.directory, self.run_id+'.yml')
        pars = self.samples[sample_number]
        coronagraph, input_pars = sli.import_pars(file_path_yml)
        (merit_val, tau, whole_throughput, throughputs, whole_contrast
         , band_contrasts, encircle_psf, wavefronts, Cr_p, Cr_b_star, Cr_b_zodi
         , Cr_b_exozodi, Cr_b) = slr.ln_probability(design_parameters=pars
                                     , coronagraph=coronagraph
                                     , non_negative=input_pars['non_negative']
                                     , input_pars=input_pars, reconstruct=True)
        return coronagraph

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
    # Show DM1
    field = wavefronts['dm1'].electric_field
    grid_scale = coronagraph.pupil_diameter
    quad_plot(field, 'DM1', grid_scale, 'x [pupil diameter]')
    # Show ON_DM2
    field = wavefronts['on_dm2'].electric_field
    grid_scale = coronagraph.pupil_diameter
    quad_plot(field, 'ON DM2', grid_scale, 'x [pupil diameter]') 
    # Show ON_DM2
    field = wavefronts['dm2'].electric_field
    grid_scale = coronagraph.pupil_diameter
    quad_plot(field, 'DM2', grid_scale, 'x [pupil diameter]') 
    # Show RELAY_PUPIL
    field = wavefronts['relay_pupil'].electric_field
    grid_scale = coronagraph.pupil_diameter
    quad_plot(field, 'RELAY PUPIL', grid_scale, 'x [pupil diameter]') 
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
#    show_run_stats(filename='20210624t1612.hdf', plot=False, long=False)    
    reconstruct('20210624t1612', 309877)
#    reconstruct('20200508t1759', 35937)
#    reconstruct('20200516t0007', 10712)
#    reconstruct('20200516t0007', 10732)
#    reconstruct('20200516t0007', 27391)
#    reconstruct('20200518t0951', 7211)
#    reconstruct('20200706t2300', 12942)
#    reconstruct('20200726t0120', 24095)
#    reconstruct('20200726t1539', 33301)
#    reconstruct('20210111t1609', 10419)
#    reconstruct('20210111t1609', 58536)
