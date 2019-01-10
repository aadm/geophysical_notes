'''
===================
aaplot.py
===================

Functions to plot stuff.

Created May 2015 by Alessandro Amato del Monte (alessandro.adm@gmail.com)

HISTORY
2015-05-01 first public release.
'''

import numpy as np
import matplotlib.pyplot as plt
import agilegeo
from pandas import Series, DataFrame
import pandas as pd

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def ampspec(signal,sr):
    '''
    Calculates amplitude spectrum of a signal with FFT.
    '''
    from scipy.fftpack import fft, fftfreq
    from scipy.interpolate import interp1d
    SIGNAL = fft(signal)
    freq = fftfreq(signal.size, d=sr)
    # Chop off the negative frequencies
    keep = freq>=0
    SIGNAL = np.abs(SIGNAL[keep])
    freq = freq[keep]
    freq0=np.linspace(freq.min(),freq.max()/2,freq.size*10)
    f = interp1d(freq, SIGNAL, kind='cubic')
    return f(freq0),freq0

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_wavelet(wavelet,time):
    '''
    Plots wavelet.

    Required timescale can be calculated with:
        >>> time=np.arange(-duration/2, duration/2 , dt)
    where duration and dt (sample rate) are same inputs given to wavelet calculation.
    '''
    plt.figure(figsize=(8,5))
    plt.plot(time,wavelet,lw=2,color='black')
    plt.fill_between(time,wavelet,0,wavelet>0.0,interpolate=False,hold=True,color='blue', alpha = 0.5)
    plt.fill_between(time,wavelet,0,wavelet<0.0,interpolate=False,hold=True,color='red', alpha = 0.5)
    plt.grid()
    plt.xlim(-0.1,0.1)
    locs,labels = plt.xticks()
    plt.xticks(locs[:-1], map(lambda x: "%d" % x, locs[:-1]*1000))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_wavelet_spectrum(wavelet,time,dt):
    '''
    Plots wavelet and its amplitude spectrum.

    Required timescale can be calculated with:
        >>> time=np.arange(-duration/2, duration/2 , dt)
    where duration and dt (sample rate) are same inputs given to wavelet calculation.
    '''
    wavelet_fft,wavelet_freq=ampspec(wavelet,dt)

    f,ax=plt.subplots(2)
    ax[0].plot(time,wavelet,lw=2,color='black')
    ax[0].grid()
    ax[1].plot(wavelet_freq,wavelet_fft,lw=2,color='black')
    ax[1].grid()
    ax[1].set_xlim(0,250)
    # from scipy.signal import argrelmax
    # peak_freq=wavelet_freq[argrelmax(wavelet_fft)][0]
    # ax[1].set_xlim(0,peak_freq*4)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_rock_grid(data,zz=1):
    '''
    Plots rock model created with make_wedge.

    INPUT
    data: 2D numpy array containing values from 1 to 3
    zz: vertical sample rate in depth
    '''
    import matplotlib.cm as cm
    cc=cm.get_cmap('copper_r',3)
    plt.figure(figsize=(12,6))
    plt.imshow(data,extent=[0,data.shape[1],data.shape[0]*zz,0],cmap=cc,interpolation='none',aspect='auto')
    cbar=plt.colorbar()
    cbar.set_ticks(range(1,4)); cbar.set_ticklabels(range(1,4))
    plt.grid()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_density(data,zz=1,seismic=True,lag=0):
    '''
    Density plot of generic 2D numpy array (seismic or any property e.g., velocity).

    INPUT
    data: 2D numpy array containing seismic or elastic property
    zz: vertical sample rate in depth or time
    seismic: True to use red-blue colorscale
    lag: lagtime at the top of the data
    '''
    plt.figure(figsize=(12,6))
    if seismic==True:
        # clip=np.amax(abs(data))
        clip=abs(np.percentile(data, 0.999))
        plt.imshow(data,extent=[0,data.shape[1],data.shape[0]*zz+lag,0+lag],cmap='RdBu',vmax=clip,vmin=-clip,aspect='auto')
    else:
        # plt.imshow(data,extent=[0,data.shape[1],data.shape[0]*zz+lag,0+lag],cmap='PiYG',aspect='auto')
        plt.imshow(data,extent=[0,data.shape[1],data.shape[0]*zz+lag,0+lag],cmap='nipy_spectral',aspect='auto')
    plt.colorbar(), plt.grid()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_wiggle(data,zz=1,skip=1,gain=1,alpha=0.7,black=False):
    '''
    Wiggle plot of generic 2D numpy array.

    INPUT
    data: 2D numpy array
    zz: vertical sample rate in depth or time
    skip: interval to choose traces to draw
    gain: multiplier applied to each trace
    '''
    [n_samples,n_traces]=data.shape
    t=range(n_samples)
    plt.figure(figsize=(9.6,6))
    for i in range(0, n_traces,skip):
        trace=gain*data[:,i] / np.max(np.abs(data))
        plt.plot(i+trace,t,color='k', linewidth=0.5)
        if black==False:
            plt.fill_betweenx(t,trace+i,i, where=trace+i>i, facecolor=[0.6,0.6,1.0], linewidth=0)
            plt.fill_betweenx(t,trace+i,i, where=trace+i<i, facecolor=[1.0,0.7,0.7], linewidth=0)
        else:
            plt.fill_betweenx(t,trace+i,i, where=trace+i>i, facecolor='black', linewidth=0, alpha=alpha)
    locs,labels=plt.yticks()
    plt.yticks(locs,[n*zz for n in locs.tolist()])
    plt.grid()
    plt.gca().invert_yaxis()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plot_partial_stacks(near,mid,far,zz=1,label=''):
    '''
    Density plot of near, mid, far stacks.

    INPUT
    near, mid, far: 2D numpy arrays containing seismic
    zz: vertical sample rate in twt
    label
    '''
    clip=np.amax([abs(near), abs(mid), abs(far)])
    f, ax = plt.subplots(nrows=1, ncols=3, figsize=(15,5))
    im0=ax[0].imshow(near,extent=[0,near.shape[1],near.shape[0]*zz,0],cmap='RdBu',vmax=clip,vmin=-clip,aspect='auto')
    ax[0].set_title(label+' (NEAR)',fontsize='small')
    im1=ax[1].imshow(mid,extent=[0,near.shape[1],near.shape[0]*zz,0],cmap='RdBu',vmax=clip,vmin=-clip,aspect='auto')
    ax[1].set_title(label+' (MID)',fontsize='small')
    im2=ax[2].imshow(far,extent=[0,near.shape[1],near.shape[0]*zz,0],cmap='RdBu',vmax=clip,vmin=-clip,aspect='auto')
    ax[2].set_title(label+' (FAR)',fontsize='small')
    ax[0].set_ylabel('twt [s]')
    cax = f.add_axes([0.925, 0.25, 0.02, 0.5])
    cbar=f.colorbar(im0, cax=cax, orientation='vertical')
    for i in range(len(ax)):
        ax[i].grid()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def update_xlabels(min_thickness,max_thickness,n_traces):
    '''
    Updates x_labels with actual thickness of model (in meters).
    '''
    locs,labels=plt.xticks()
    incr=(max_thickness-min_thickness)/float(n_traces)
    newlabels=(locs[1:-1])*incr+min_thickness
    plt.xticks(locs[1:-1],[str(round(x,1))+'m' for x in newlabels])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def update_ylabels(lag,thickness,vel):
    '''
    Updates y_labels to add lag in two-way-time,
    given velocity of top layer having certain thickness.
    '''
    locs,labels=plt.yticks()
    lagtop=thickness/vel*2
    plt.yticks(locs[:-1],[round(y+lag-lagtop,3) for y in locs])


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def wellplot(ll,ztop=None,zbot=None):
    '''
    wellplot (c) aadm 2014
    Simple logview plot.
    Needs in input a Pandas dataframe containing log structure, and optionally depth range.

    HISTORY
    2014-08-04 updated with depth range input
    2014-06-30 first version
    '''
    if ztop==None:
        ztop=ll.DEPTH.min()
    if zbot==None:
        zbot=ll.DEPTH.max()
    velmin=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['VS']].min().values
    velmax=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['VP']].max().values
    dmin=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['RHO']].min().values
    dmax=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['RHO']].max().values
    ipmin=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['IP']].min().values
    ipmax=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['IP']].max().values
    rmin=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['VPVS']].min().values
    rmax=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['VPVS']].max().values

    swplot=True if 'SW' in ll.columns else False
    f, ax = plt.subplots(nrows=1, ncols=5, sharey=True, figsize=(12,6))
    if swplot:
        ll.plot(x='SW', y='DEPTH',ax=ax[0], style='b', label='Sw');
    ll.plot(x='VSH',  y='DEPTH',  ax=ax[0], style='g', label='Vsh');
    ll.plot(x='PHI',  y='DEPTH',  ax=ax[0], style='k', label='phi');
    ll.plot(x='VP',   y='DEPTH',  ax=ax[1], style='k');
    ll.plot(x='VS',   y='DEPTH',  ax=ax[1], style='r');
    ll.plot(x='RHO',  y='DEPTH',  ax=ax[2], style='k');
    ll.plot(x='IP',   y='DEPTH',  ax=ax[3], style='k');
    ll.plot(x='VPVS', y='DEPTH',  ax=ax[4], style='k');
    ax[0].set_xlabel('Vcl/Sw/phi'),       ax[0].set_xlim(-0.1,1.1),                          ax[0].set_ylim(ztop,zbot)
    ax[1].set_xlabel('Velocities (m/s)'), ax[1].set_xlim(velmin-velmin*.1,velmax+velmax*.1), ax[1].set_ylim(ztop,zbot)
    ax[2].set_xlabel('Density (g/cc)'),   ax[2].set_xlim(dmin-dmin*.01,dmax+dmax*.01),       ax[2].set_ylim(ztop,zbot)
    ax[3].set_xlabel('Ip (m/s*g/cc)'),    ax[3].set_xlim(ipmin-ipmin*.1,ipmax+ipmax*.1),     ax[3].set_ylim(ztop,zbot)
    ax[4].set_xlabel('Vp/Vs'),            ax[4].set_xlim(rmin-rmin*.01,rmax+rmax*.01),       ax[4].set_ylim(ztop,zbot)
    ax[0].invert_yaxis()
    for i in range(len(ax)):
        ax[i].locator_params(axis='x', nbins=4)
    ax[0].legend(fontsize='small', loc='lower left')
    ax[1].legend(fontsize='small', loc='lower left')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def plotlfc(x,y,z,bb=30,colormap='Paired'):
    '''
    plotlfc (c) aadm 2014
    Crossplot of x vs y color coded to z.
    Useful for crossplotting two log properties color coded to a LFC (litho-fluid class, or facies).
    Plots the same data twice, first as scatterpoints then as contours after calculating a 2D histogram.

    Optional parameter bb controls 2D histogram bins (default=30).

    HISTORY
    2014-05-19 first version
    '''
    import matplotlib.cm
    import matplotlib.colors

    nfacies=len(np.unique(z))
    #.......................... qui di seguito un accrocchio per discretizzare la colormap
    colori = plt.get_cmap(colormap, nfacies)
    cNorm=matplotlib.colors.Normalize(vmin=1, vmax=nfacies+1)
    scalarMap=matplotlib.cm.ScalarMappable(norm=cNorm,cmap=colori)
    ccc=scalarMap.to_rgba(range(1,nfacies+1))
    #......................................................................................
    fig, ax = plt.subplots(1, figsize=(8,6))
    ax.scatter(x, y, 20, z, cmap=colori, alpha=0.7, marker='o', edgecolors='none')
    colorbar_index(ncolors=nfacies, cmap=colori)
    dimfig=[ax.get_position().x0, ax.get_position().y0, ax.get_position().x1, ax.get_position().y1]

    fig, ax = plt.subplots(1, figsize=(8,6))
    for i in np.unique(z):
        H, xedges, yedges = np.histogram2d(x[z==i], y[z==i], bins=bb) # normed=False (default), returns the number of samples in each bin.
        H = blur_image(H,3)
        xi = np.linspace(np.min(x[z==i]), np.max(x[z==i]), H.shape[0])
        yi = np.linspace(np.min(y[z==i]), np.max(y[z==i]), H.shape[0])
        cset=plt.contour(xi, yi, H, 2, linewidths=2, colors=matplotlib.colors.rgb2hex(ccc[int(i)-1]), label=('LFC=%d' %i))
    plt.subplots_adjust(dimfig[0],dimfig[1],dimfig[2],dimfig[3])

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def colorbar_index(ncolors, cmap):
    '''
    Found somewhere on stackoverflow,
    useful to get a colorbar for discrete colors with values written in the middle of each color patch.
    '''
    mappable = plt.cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolors+0.5)
    colorbar = plt.colorbar(mappable)
    colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
    colorbar.set_ticklabels(range(0,ncolors+1))
