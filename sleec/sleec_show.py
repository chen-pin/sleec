import os
import pickle as pkl
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
    iteration : int
        Number of iterations

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

    def summary(self):
        """ 
        Print summary of EMCEE run characteristics
        """
        log_prob_max = self.log_prob_samples.max()
        ndx = np.argmax(self.log_prob_samples)
        print("number of walkers: {0}".format(self.num_walkers))
        print("number of dimensions: {0}".format(self.num_dims))
        print("number of iterations: {0}".format(self.iteration))
        print("acceptance fraction:  {0}".format(self.acceptance_fractn()))
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

    def top_samples(self, n=5, print_summary=True, print_pars=False):
        """
        Return `self.samples` indices of the `n` samples with the highest 
        ln probability values.  

        """
        sortarg = np.flip(np.argsort(self.log_prob_samples))
        ndx = [sortarg[0]]
        count = 1 
        prev_log_prob = self.log_prob_samples[sortarg[0]]
        if print_summary:
            for i in sortarg[1:]:
                log_prob = self.log_prob_samples[i]
                if (log_prob != prev_log_prob and count<=n-1):
                    ndx.append(i)
                    prev_log_prob = log_prob
                    count += 1
#                elif count>3:
#                    break
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
            for j in ndx:
               print("\nparameter values:  {}\n"
                         .format(self.samples[j]))
        return ndx[:n]

    def corner_plot(self, discard=0, save=False):
        all_samples = np.concatenate(
                (self.samples[discard:], self.log_prob_samples[discard:, None])
                 , axis=1)
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

    def reconstruct(self, sample_number, save, verbose=True, show=True):  
        file_path_slp = os.path.join(self.directory, self.run_id+'.pkl')
        file_path_yml = os.path.join(self.directory, self.run_id+'.yml')
        pars = self.samples[sample_number]
        coronagraph, input_pars = sli.import_pars(file_path_yml)
        (coronagraph, ln_merit, tau, whole_throughput, throughputs
         , whole_contrast, band_contrasts, encircle_psf, wavefronts, Cr_p
         , Cr_b_star, Cr_b_zodi, Cr_b_exozodi
         , Cr_b) = slr.ln_probability(design_parameters=pars
                                     , coronagraph=coronagraph
                                     , non_negative=input_pars['non_negative']
                                     , input_pars=input_pars, reconstruct=True)
        if verbose:
            print("ln merit value:  {}".format(ln_merit))
            print("required exposure time:  {}".format(tau))
            print("overall throughput:  {}".format(whole_throughput))
            print("throughput each Airy band:  {}".format(throughputs))
            print("overall contrast:  {}".format(whole_contrast))
            print("contrast in each Airy band:  {}".format(band_contrasts))
            print("encircled PSF energy {}".format(encircle_psf))
        if (show or save):
            iwa = input_pars['IWA']
            owa = input_pars['OWA']
            resel = coronagraph.f_num*coronagraph.wvl
            psf = wavefronts['on_fpm'].intensity
            psf0 = psf.max()
            radii = wavefronts['on_fpm'].electric_field.grid.as_('polar').r
            mask = np.where((iwa*resel <= radii) & (radii < owa*resel), 1, 0)
            avg_background = (wavefronts['camera'].intensity[np.where(mask)]
                                                       .mean())
            contrast_r = avg_background/psf0
            print("recomputed contrast:  {:E}".format(contrast_r))
            # Show PUPIL
            field = coronagraph.pupil_field
            field = field/np.abs(field).max()
            grid_scale = coronagraph.pupil_diameter/2.0
            if save:
                fn = 'RUN'+self.run_id+'_SAMPLE'+str(sample_number)+'_PUPIL.png'
                path = os.path.join(self.directory, fn)
                quad_plot(field, 'PUPIL', grid_scale, 'x [pupil radius]'
                          , show=show, scale1='linear', save=save, path=path)
            else:
                quad_plot(field, 'PUPIL', grid_scale, 'x [pupil radius]'
                          , show=show, scale1='linear', save=save)
            
            # Show DM1_SURFC
            field = coronagraph.dm1_surfc
            if save:
                fn = 'RUN'+self.run_id+'_SAMPLE'+str(sample_number)+'_DM1-SURFC.png'
                path = os.path.join(self.directory, fn)
                single_plot(field, title='DM1 SURFACE'
                            , xlabel ='x [m]', scale='linear'
                            , vmin=field.min(), vmax=field.max()
                            , show=show, save=save, path=path)
            else:
                single_plot(field, title='DM1 SURFACE'
                            , xlabel ='x [m]', scale='linear'
                            , vmin=field.min(), vmax=field.max()
                            , show=show, save=save)
            # Show DM1
            field = wavefronts['dm1'].electric_field
            if save:
                fn = 'RUN'+self.run_id+'_SAMPLE'+str(sample_number)+'_POST-DM1.png'
                path = os.path.join(self.directory, fn)
                quad_plot(field, 'POST-DM1', grid_scale, 'x [pupil radius]'
                          , show=show, scale1='linear', save=save, path=path)
            else:
                quad_plot(field, 'POST-DM1', grid_scale, 'x [pupil radius]'
                          , show=show, scale1='linear', save=save)
            # Show ON_DM2
            field = wavefronts['on_dm2'].electric_field
            field = field/np.abs(field).max()
            if save:
                fn = 'RUN'+self.run_id+'_SAMPLE'+str(sample_number)+'_ON-DM2.png'
                path = os.path.join(self.directory, fn)
                quad_plot(field, 'ON DM2', grid_scale, 'x [pupil radius]'
                          , show=show, scale1='linear', save=save, path=path)
            else:
                quad_plot(field, 'ON DM2', grid_scale, 'x [pupil radius]'
                          , show=show, scale1='linear', save=save)
            # Show DM2_SURFC
            field = coronagraph.dm2_surfc
            if save:
                fn = 'RUN'+self.run_id+'_SAMPLE'+str(sample_number)+'_DM2-SURFC.png'
                path = os.path.join(self.directory, fn)
                single_plot(field, title='DM2 SURFACE'
                            , xlabel ='x [m]', scale='linear'
                            , vmin=field.min(), vmax=field.max()
                            , show=show, save=save, path=path)
            else:
                single_plot(field, title='DM2 SURFACE'
                            , xlabel ='x [m]', scale='linear'
                            , vmin=field.min(), vmax=field.max()
                            , show=show, save=save)
            # Show DM2
            field = wavefronts['dm2'].electric_field
            if save:
                fn = 'RUN'+self.run_id+'_SAMPLE'+str(sample_number)+'_POST-DM2.png'
                path = os.path.join(self.directory, fn)
                quad_plot(field, 'POST-DM2', grid_scale, 'x [pupil radius]'
                          , show=show, scale1='linear', save=save, path=path)
            else:
                quad_plot(field, 'POST-DM2', grid_scale, 'x [pupil radius]'
                          , show=show, scale1='linear', save=save)
            # Show RELAY_PUPIL
            field = wavefronts['relay_pupil'].electric_field
            field = field/np.abs(field).max()
            if save:
                fn = 'RUN'+self.run_id+'_SAMPLE'+str(sample_number)+'_RELAY-PUPIL.png'
                path = os.path.join(self.directory, fn)
                quad_plot(field, 'RELAY PUPIL', grid_scale, 'x [pupil radius]'
                          , show=show, scale1='linear', save=save, path=path)
            else:
                quad_plot(field, 'RELAY PUPIL', grid_scale, 'x [pupil radius]'
                          , show=show, scale1='linear', save=save)
            # Show ON_FPM
            field = wavefronts['on_fpm'].electric_field
            field = field/np.abs(field).max()
            grid_scale = coronagraph.f_num*coronagraph.wvl
            scale_max = np.ceil(np.abs(field.max())**2) 
            scale_min = 1e-13
            if save:
                fn = 'RUN'+self.run_id+'_SAMPLE'+str(sample_number)+'_ON-FPM.png'
                path = os.path.join(self.directory, fn)
                quad_plot(field, 'INCIDENT FIELD ON THE FPM', grid_scale
                          , '[lambda/D]', psf_vmin=scale_min, psf_vmax=scale_max
                          , show=show, scale1='log', scale2='log'
                          , save=save, path=path)
            else:
                quad_plot(field, 'INCIDENT FIELD ON THE FPM', grid_scale
                          , '[lambda/D]', psf_vmin=scale_min, psf_vmax=scale_max
                          , show=show, scale1='log', scale2='log'
                          , save=save )
            # Show FPM
            field = coronagraph.fpm_field
            grid_scale = coronagraph.f_num*coronagraph.wvl
            if save:
                fn = 'RUN'+self.run_id+'_SAMPLE'+str(sample_number)+'_FPM.png'
                path = os.path.join(self.directory, fn)
                quad_plot(field, 'FPM', grid_scale, '[lambda/D]'
                          , show=show, scale1='linear', save=save, path=path)
            else:
                quad_plot(field, 'FPM', grid_scale, '[lambda/D]'
                          , show=show, scale1='linear', save=save)
            # Show ON_LYOT
            field = wavefronts['on_lyot'].electric_field
            field = field/np.abs(field).max()
            if save:
                fn = 'RUN'+self.run_id+'_SAMPLE'+str(sample_number)+'_ON-LYOT.png'
                path = os.path.join(self.directory, fn)
                quad_plot(field, 'INCIDENT FIELD ON THE LYOT STOP', grid_scale, 'x [pupil radius'
                          , show=show, scale1='linear', save=save, path=path)
            else:
                quad_plot(field, 'INCIDENT FIELD ON THE LYOT STOP', grid_scale, 'x [pupil radius]'
                          , show=show, scale1='linear', save=save)
            # Show LYOT APERTURE
            field = coronagraph.lyot_field
            if save:
                fn = 'RUN'+self.run_id+'_SAMPLE'+str(sample_number)+'_LYOT.png'
                path = os.path.join(self.directory, fn)
                quad_plot(field, 'LYOT STOP', grid_scale, 'x [pupil radius]'
                          , show=show, scale1='linear', save=save, path=path)
            else:
                quad_plot(field, 'LYOT STOP', grid_scale, 'x [pupil radius]'
                          , show=show, scale1='linear', save=save)
            # Show CONTRAST	
            field = ( np.abs(wavefronts['camera'].electric_field)
                      /np.abs(wavefronts['on_fpm'].electric_field).max() )
            scale_max = 1e-9
            scale_min = 1e-13
            grid_scale = coronagraph.f_num*coronagraph.wvl
            if save:
                fn = 'RUN'+self.run_id+'_SAMPLE'+str(sample_number)+'_CONTRAST.png'
                path = os.path.join(self.directory, fn)
                quad_plot(field, 'CONTRAST', grid_scale, '[lambda/D]'
                          , normalization='none'
                          , scale='log', psf_vmin=scale_min, psf_vmax=scale_max
                          , show=show, save=save, path=path) 
            else:
                quad_plot(field, 'CONTRAST', grid_scale, '[lambda/D]'
                          , scale='log', psf_vmin=scale_min, psf_vmax=scale_max
                          , normalization='none'
                          , show=show) 
        return_obj = (coronagraph, ln_merit, tau, whole_throughput, throughputs
            , whole_contrast, band_contrasts, encircle_psf, wavefronts, Cr_p
            , Cr_b_star, Cr_b_zodi, Cr_b_exozodi, Cr_b)  
#        if save:
#            fn = 'RUN'+self.run_id+'_SAMPLE'+str(sample_number)+'_RETURN.pkl'
#            path = os.path.join(self.directory, fn)
#            with open(path, 'wb') as f:
#                pkl.dump(return_obj, f)

        return return_obj
        
    def acceptance_fractn(self):
        return np.mean(self.accepted/float(self.iteration))

def quad_plot(field, suptitle, grid_scale, xlabel, scale1='linear'
              , scale2='linear', scale3='linear', scale4='linear'
#              , normalization='peak'
              , show=True, psf_vmin=1E-12, psf_vmax=1.0
              , save=False, **kwargs):
    intensity = (np.abs(field))**2
    if np.iscomplexobj(field):
        phase = np.angle(field)
    else:
        phase = np.zeros_like(field)
    plt.figure(figsize=(15, 12))
    plt.suptitle(suptitle)
    plt.subplot(221)
    plt.title('intensity')
    hci.imshow_psf(intensity, scale=scale1
#                   , normalization=normalization
                   , cmap=turbo, vmin=psf_vmin, vmax=psf_vmax
                   , spatial_resolution=grid_scale)
    plt.subplot(223)
    plt.title('phase [pi]')
    hci.imshow_psf(phase/np.pi, scale=scale3
                   , cmap=turbo, vmax=1.0, vmin=-1.0
                   , spatial_resolution=grid_scale)
    plt.subplot(222)
    side_dim = int(np.sqrt(len(field)))
    ctr_ndx = side_dim//2
    plt.title('Horizontal Center Cross-Section:  Intensity')
    abscissa = ( field.grid.x.reshape(side_dim, side_dim)
                                                     [ctr_ndx, :] )/grid_scale
    plt.yscale(scale2)
    plt.plot( abscissa 
              , intensity.reshape(side_dim, side_dim)[ctr_ndx,:])
    plt.xlabel(xlabel)
#    plt.ylabel('amplitude')
    plt.subplot(224)
    plt.title('Horizontal Center Cross-Section:  Phase')
    plt.yscale(scale4)
    plt.plot( abscissa 
              , phase.reshape(side_dim, side_dim)[ctr_ndx,:]/np.pi)
    plt.xlabel(xlabel)
    plt.ylabel('phase [pi]')
    if save:
        plt.savefig(kwargs['path'])
    if show:
        plt.show()

def single_plot(field, title, xlabel, scale='linear'
                , vmin=1e-12, vmax=1.0
                , show=True, save=False, **kwargs):
    plt.figure(figsize=(15, 12))
    plt.title(title)
    hci.imshow_psf(field, scale=scale, cmap=turbo, vmin=vmin, vmax=vmax)
    plt.xlabel(xlabel)
    if save:
        plt.savefig(kwargs['path'])
    if show:
        plt.show()


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
