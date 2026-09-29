def tilted_circular_aperture(diameter, x_tilt=0.0, y_tilt=0.0, x_shift=0.0
                             , y_shift=0.0):
    """Makes a Field generator for a circular aperture with linear phase shift.

    Parameters
    ----------
    diameter : scalar
        The diameter of the aperture.
    tilt : scalar 
        Phase shift (radians) across the diameter.
    center : array_like
        The center of the aperture
     
    Returns
    -------
    Field generator
        This function can be evaluated on a grid to get a Field
   
    Notes
    -----
    This function is Coded but not testd, likely has bugs.  It is a revision of 
    John Trauger's function by the same name.  
   
    """

    def func(grid):
        polar =  grid.as_('polar')
        f = np.zeros_like(polar.r,'complex128')
        f[np.where(polar.r <= diameter/2)] = 1.0
        f *= np.exp( (1j*x_tilt*(x-x_shift) + 1j*y_tilt*(y-y_shift) 
                     / diameter)
            else:
                f = grid.r <= (diameter / 2)
        return Field(f.astype('complex128'), grid)
    return func


def update_pupil(coro, y_tilt_airy, x_tilt_airy):
    ytilt = 2.0*np.pi*y_tilt_airy   # tilt in radians 
    xtilt = 2.0*np.pi*x_tilt_airy   # tilt in radians 
    pupil = hci.tilted_circular_aperture(diameter=coro.pupil_diamtr
                                         , ytilt=ytilt, xtilt=xtilt)
                                        (coro.pupil_field.grid)
    coro.pupil_field = pupil
    coro.optical_elements[0] = coro.pupil
    return coro


def planetPSF(coro, incdnt_wavfrnt, starMax, sumStarPSF, tilt_airy):
    # planet PSF returns number of pixels within the FWHM  and pixel 
    # intensities summed 
    #     over the FWHM, for the offset star in the coro (camera) focal plane
    update_pupil(coro, y_tilt_airy=tilt_airy, x_tilt_airy=0.0)
    output_tilted = coro(incdnt_wavfrnt)
    PSF = output_tilted['camera'].intensity
    mask = np.zeros_like(PSF)
    ii = np.where(PSF >= 0.5*PSF.max())
    mask[ii] = 1
    sumPSF = np.sum(mask*PSF)
    pixPSF = np.sum(mask)
    hwhmPSF = np.sqrt(pixPSF/np.pi) / coro.pix_per_airy  # define annular 
        # scoring footprint
    r1 = (tilt_airy - hwhmPSF)*coro.f_num*coro.wvl
    r2 = (tilt_airy + hwhmPSF)*coro.f_num*coro.wvl
    update_pupil(coro, y_tilt_airy=0.0, x_tilt_airy=0.0)
    output_untilted = coro(incdnt_wavfrnt)
    pixDonut, sumDonut = sum_dark_annulus(r1, r2, coro, output_untilted)
    makefits = False
    if makefits:
        imgdim = int(np.sqrt(PSF.shape))
        img = np.reshape(PSF,[imgdim,imgdim])
        fits3.wfits("planet_image.fits",img)
#        fits3.wfits("donut_image.fits",img+dmask*img)
    return  pixPSF, sumPSF, pixDonut, sumDonut

