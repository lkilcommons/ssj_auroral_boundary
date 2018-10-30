"""
Some simple code to make particle flux spectrograms with matplotlib
@author: Liam M. Kilcommons
		 (minor modifications R. Redmon)
"""
import numpy as np
import matplotlib.pyplot as pp
import datetime as dt

def dmsp_spectrogram(times, flux, channel_energies=None, lat=None, lt=None,
                     fluxunits='eV/cm$^2$-s-sr-eV', logy=True, datalabel=None,
                     cblims=None, title=None, ax=None, ax_cb=None,
                     label_it=True):
    """ Plot the DMSP spectrogram

    Parameters
    ----------
    times : numpy.ndarray (dtype=object)(shape=(n,1))
        Array of datetimes corresponding to the timestamps of the rows of the
        flux array
    flux : numpy.ndarray (shape=(n,len(channel_energies)))
        Array of fluxes, 1 per channel, per timestamp
    channel_energies - numpy.ndarray
        Array of particle detector channel center energies in eV, if 
        None uses default DMSP energies
        channel_energies = [ 30000.,  20400.,  13900.,   9450.,   6460.,
                              4400.,   3000.,   2040.,   1392.,    949.,
                               646.,    440.,    300.,    204.,    139.,
                                95.,     65.,     44.,     30.]
    fluxunits : str, optional
        Units of flux for labeling the spectrogram (title and colorbar)
        Defaults to eV/cm$^2$-s-sr-eV
    logy : boolean, optional
        Flag to make the y axis log scale
        (useful for log-spaced channel_energies)
    lat : numpy.ndarray (shape=(n,1)), optional
        If lat is not None, then it must be the latitude
        (magnetic or otherwise) of the spacecraft at
        every timestamp in times. Setting this value
        will cause the latitude to be added to the
        x axis labels
    lt : numpy.ndarray (shape=(n,1)), optional
        If lat is not None, then it must be the localtime
        (magnetic or otherwise) of the spacecraft at
        every timestamp in times. Setting this value
        will cause the localtime to be added to the
        x axis labels
    datalabel : str, optional
        Some text to add to the title of the graphic
        goes on a line above 'Flux [units of flux]'
    cblims : None or 2-element list, optional
        The limits for the colorbar. If None,
        then the colorbar range is set to [flux.min(),flux.max()]
    ax : None or axis reference, optional
        Allows caller to specify axis for spectrogram; helpful for stackplot.
        If 'ax' is specified then so should 'ax_cb'.
    ax_cb : None or colorbar axis reference, optional
        Allows caller to specify axis for spectrogram color bar; helpful for
        stackplot. If 'ax' is specified then so should 'ax_cb'.

    """
    #Module for logrithmic colorbar spacing
    from matplotlib.colors import LogNorm
    #Module for locating dates on the x axis
    import matplotlib.dates as mpldates
    #Module for getting colormaps
    import matplotlib.cm as cm

    if channel_energies is None:
        channel_energies = np.array([ 30000., 20400., 13900., 9450., 6460.,
                                      4400., 3000., 2040., 1392., 949., 646.,
                                      440., 300., 204., 139., 95., 65., 44.,
                                      30.])

    # if Axis not specified then create one
    if ax is None:
        f = pp.figure(figsize=(12,6),dpi=300)
        ax = pp.axes()

    if datalabel is not None:
        ax.set_title(datalabel+'\n Flux [%s]' %(fluxunits))
    else:
        pass
        #ax.set_title('Flux [%s]' % (fluxunits))
    if isinstance(times,np.ndarray):    
        times = times.flatten()
    
    if isinstance(times[0], dt.datetime):
        mpl_times = mpldates.date2num(times)
    else:
        mpl_times = times

    #--------------------------------------------------------------------------
    # Channel center energies to bin starts
    # Since DMSP SSJ channels are log-linearly spaced, the bins widths are taken
    # to be log-constant and the bins are placed symmetric about the channel
    # center energies. This is probably not exactly correct since the actual
    # instrument response/sensitivity functions are likely more linear than
    # log linear. Recall that channels are listed as [30,000 eV to 30 eV] in
    # reverse order.
    #--------------------------------------------------------------------------

    # Hard coded start/end bins taken from SSDP; not sure how they are derived,
    # though this does result in bins visually centered correctly on their
    # central energies
    bin_edges = np.logspace(np.log10(36340.), np.log10(24.76),
                            len(channel_energies) + 1) # add one for endpoint
    T,CH_E = np.meshgrid(mpl_times, bin_edges)

    # Infinite, and Negative fluxes => NaN
    inds = np.nonzero((~np.isfinite(flux)) | (flux < 0.))
    flux[inds] = np.nan

    # Mask nan fluxes so that pcolor knows to use the cmap bad value
    masked_flux = np.ma.masked_where(np.isnan(flux),flux)

    if cblims is None:
        z_min = np.nanmin(flux)
        z_max = np.nanmax(flux)
    else:
        z_min = cblims[0]
        z_max = cblims[1]

    #Set the over and under-range colors for the colorbar
    cmap = cm.get_cmap('jet')
    cmap.set_bad('white',.1)
    cmap.set_over('black')
    cmap.set_under('grey')

    mappable = ax.pcolormesh(T, CH_E, masked_flux.transpose(), cmap=cmap,
                             norm=LogNorm(vmin=z_min, vmax=z_max))
    #mappable.set_rasterized( True )

    if ax_cb is None:
        pp.colorbar(mappable,label=fluxunits,ax=ax)
    else:
        pp.colorbar(mappable,label=fluxunits,cax=ax_cb)
        
    # if Axis not specified then add x-axis tick marks
    if label_it and isinstance(times[0], dt.datetime):

        plotwidth_h = (times[-1]-times[0]).total_seconds()/3600.
        plotwidth_m = (times[-1]-times[0]).total_seconds()/60.

        if  plotwidth_m <= 10.:
            # if the plot width is less than 10 minutes tick mark every minute
            majloc = mpldates.MinuteLocator(interval=1)
        elif  plotwidth_m <= 30.:
            # if the plot width is less than 1/2 hour tick mark every 5 minutes
            majloc = mpldates.MinuteLocator(interval=5)
        elif  plotwidth_h <= 1:
            # if the plot width is less than 1 hour, but more than 30 minutes,
            # tick mark every 10 minutes
            majloc = mpldates.MinuteLocator(interval=10)
        elif plotwidth_h <= 3:
            # if less than 3 hours, but more than 1 use every 15 minutes
            majloc = mpldates.MinuteLocator(interval=15)
        elif plotwidth_h <= 5:
            # if less than 5 hours, but more than 3 use every half hour
            majloc = mpldates.MinuteLocator(interval=30)
        else:
            majloc = mpldates.HourLocator() #tick mark every hour

        #Set the date locator
        ax.xaxis.set_major_locator(majloc)

        #This is slow and throws errors if used with pcolor, used pcolormesh
        # instead
        #ax.set_yscale('log')

        #Manually create the tick labels
        #There is probably a better way to do this with FuncFormatter, but I
        # couldn't figure out how to get all of the relavent lat and LT
        # information into it

        #Get the tick marks
        xticks = ax.get_xticks()

        xlabels = []
        for tick in xticks:
            ind = np.nonzero(mpl_times==tick)[0] #Nonzero returns array ARG!
            if len(ind)>0:
                #Sometimes tick is not found if it wants to tickmark outside of
                # data range.  Have to put additional index to get datetime
                # instead of array of length 1 with datetime in it
                tickstr = "%.2d:%.2d" % (times[ind[0]].hour,times[ind[0]].minute)
                if lat is not None:
                    tickstr+="\n%.2f" % (lat[ind])
                if lt is not None:
                    tickstr+="\n%.2f" % (lt[ind])
                xlabels.append(tickstr)
            else:
                dtime = mpldates.num2date(tick) #Convert the tick position to a time
                xlabels.append('%.2d:%.2d' % (dtime.hour, dtime.minute))

        ax.set_xticklabels(xlabels)

    ax.set_yscale('log')
    ax.set_ylim([channel_energies.min(),channel_energies.max()])
    ax.set_ylabel('Channel E \n(log)[eV]')

    # In the case that caller didn't specify the axis to use return new figure
    if 'f' in locals():
        # f.savefig('/home/liamk/test.png',dpi=300,figsize=(12,6))
        return f
