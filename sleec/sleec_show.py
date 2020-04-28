import os
import numpy as np
import matplotlib.pyplot as plt
import corner
import emcee as emc


def show_run_stats():
    filename = input('Enter EMCEE backend filename:  ')
    path = os.path.join('..', '..', 'sleec_out', filename)
    reader = emc.backends.HDFBackend(path)
    tau = reader.get_autocorr_time(tol=4)
    burnin = int(2 * np.max(tau))
    thin = int(0.5 * np.min(tau))
    samples = reader.get_chain(discard=burnin, flat=True, thin=thin)
#    samples = reader.get_chain(discard=burnin, flat=True)
    ndim = samples.shape[-1]
    log_prob_samples = reader.get_log_prob(discard=burnin, flat=True, thin=thin)
#    log_prob_samples = reader.get_log_prob(discard=burnin, flat=True)
    log_prob_max = log_prob_samples.max()
    log_prob_max_ndx = np.argmax(log_prob_samples)
    print("burn-in: {0}".format(burnin))
    print("thin: {0}".format(thin))
    print("flat chain shape: {0}".format(samples.shape))
    print("flat log prob shape: {0}".format(log_prob_samples.shape))
    print("maximum log prob value:  {}".format(log_prob_max))
    print("index of max value:  {}".format(log_prob_max_ndx))
    for i in range(ndim): 
        print("optimal parameter {} value:  {}"
              .format(i, samples[log_prob_max_ndx, i]))
    all_samples = np.concatenate(
        (samples, log_prob_samples[:, None]), axis=1
                                )
    labels = list(map(r"$\theta_{{{0}}}$".format, range(1, ndim + 1)))
    labels += ["log prob"]
    corner.corner(all_samples, labels=labels)
    plt.show()


if __name__ == "__main__":
    show_run_stats()    
