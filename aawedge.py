'''
===================
aawedge.py
===================

Functions to build and plot seismic wedges.

Created April 2015 by Alessandro Amato del Monte (alessandro.adm@gmail.com)

Heavily inspired by Matt Hall and Evan Bianco's blog posts and code:

http://nbviewer.ipython.org/github/agile-geoscience/notebooks/blob/master/To_make_a_wedge.ipynb
http://nbviewer.ipython.org/github/kwinkunks/notebooks/blob/master/Spectral_wedge.ipynb
http://nbviewer.ipython.org/github/kwinkunks/notebooks/blob/master/Faster_wedges.ipynb
http://nbviewer.ipython.org/github/kwinkunks/notebooks/blob/master/Variable_wedge.ipynb

Also see Wes Hamlyn's tutorial on Leading Edge "Thin Beds, tuning and AVO" (December 2014):

https://github.com/seg/tutorials/tree/master/1412_Tuning_and_AVO

HISTORY
2015-05-07 updated make_synth, now works also on 1D arrays.
2015-04-10 first public release.
'''

import numpy as np
import matplotlib.pyplot as plt
import agilegeo

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def make_wedge(n_traces,encasing_thickness,min_thickness,max_thickness,dz=0.1):
    '''
    Creates wedge-shaped model made of 3 units with variable thickness.

    INPUT
    n_traces
    encasing_thickness
    min_thickness
    max_thickness
    dz: vertical sample rate, by default 0.1 m

    OUTPUT
    wedge: 2D numpy array containing wedge-shaped model made of 3 units
    '''
    encasing_thickness *= (1./dz)
    min_thickness *= (1./dz)
    max_thickness *= (1./dz)
    deltaz=float(max_thickness-min_thickness)/float(n_traces)
    n_samples=max_thickness+encasing_thickness*2
    top_wedge=encasing_thickness
    wedge = np.zeros((n_samples, n_traces))
    wedge[0:encasing_thickness,:]=1
    wedge[encasing_thickness:,:]=3
    wedge[encasing_thickness:encasing_thickness+min_thickness,:]=2
    for i in range(n_traces):
        wedge[encasing_thickness+min_thickness:encasing_thickness+min_thickness+int(round(deltaz*i)),i]=2
    print "wedge minimum thickness: %.2f m" % (min_thickness*dz)
    print "wedge maximum thickness: %.2f m" % (max_thickness*dz)
    print "wedge vertical sampling: %.2f m" % (dz)
    print "wedge samples, traces: %dx%d" % (wedge.shape)
    return wedge

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def assign_ai(model, aiprop):
    '''
    Assigns acoustic impedance to a rock model created with make_wedge.

    INPUT
    model: 2D numpy array containing values from 1 to 3
    aiprop: np.array([[vp1,rho1],[vp2,rho2],[vp3,rho3]])

    OUTPUT
    model_ai: 2D numpy array containing acoustic impedances
    '''
    model_ai=np.zeros(model.shape)
    code = 1
    for x in aiprop:
        model_ai[model==code] = x[0]*x[1]
        code += 1
    return model_ai

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def assign_vel(model, aiprop):
    '''
    Assigns velocity to a rock model created with make_wedge,
    to be used for depth-time conversion.

    INPUT
    model: 2D numpy array containing values from 1 to 3
    aiprop: np.array([[vp1,rho1],[vp2,rho2],[vp3,rho3]])

    OUTPUT
    model_vel: 2D numpy array containing velocities
    '''
    model_vel=np.zeros(model.shape)
    code=1
    for x in aiprop:
        model_vel[model==code] = x[0]
        code += 1
    return model_vel

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def assign_el(model, elprop):
    '''
    Assigns elastic properties (Vp, Vs, rho) to a rock model created with make_wedge.

    INPUT
    model: 2D numpy array containing values from 1 to 3
    elprop: np.array([[vp1,rho1,vs1],[vp2,rho2,vs2],[vp3,rho3,vs3]])

    OUTPUT
    model_vp: 2D numpy array containing Vp
    model_vs: 2D numpy array containing Vs
    model_rho: 2D numpy array containing densities
    '''
    model_vp=np.zeros(model.shape)
    model_vs=np.zeros(model.shape)
    model_rho=np.zeros(model.shape)
    code = 1
    for i in elprop:
        model_vp[model==code]  = i[0]
        model_vs[model==code]  = i[2]
        model_rho[model==code] = i[1]
        code += 1
    return model_vp,model_vs,model_rho

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def make_rc(model_ai):
    '''
    Computes reflectivities of an acoustic model created with make_wedge + assign_ai.

    INPUT
    model: 2D numpy array containing acoustic impedances

    OUTPUT
    rc: 2D numpy array containing reflectivities
    '''
    upper = model_ai[:-1][:][:]
    lower = model_ai[1:][:][:]
    rc=(lower - upper) / (lower + upper)
    if model_ai.ndim==1:
        rc=np.concatenate((rc,[0]))
    else:
        n_traces=model_ai.shape[1]
        rc=np.concatenate((rc,np.zeros((1,n_traces))))  # add 1 row of zeros at the end
    return rc

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def make_rc_elastic(model_vp,model_vs,model_rho,ang):
    '''
    Computes angle-dependent reflectivities of an elastic model created with make_wedge + assign_el.
    Uses Aki-Richards approximation.

    INPUT
    model_vp: 2D numpy array containing Vp values
    model_vs: 2D numpy array containing Vs values
    model_rho: 2D numpy array containing density values
    ang: list with near, mid, far angle, e.g. ang=[5,20,40]

    OUTPUT
    rc_near: 2D numpy array containing near-stack reflectivities
    rc_mid: 2D numpy array containing mid-stack reflectivities
    rc_far: 2D numpy array containing far-stack reflectivities
    '''
    from agilegeo.avo import akirichards
    [n_samples, n_traces] = model_vp.shape
    rc_near=np.zeros((n_samples,n_traces))
    rc_mid=np.zeros((n_samples,n_traces))
    rc_far=np.zeros((n_samples,n_traces))
    uvp  = model_vp[:-1][:][:]
    lvp  = model_vp[1:][:][:]
    uvs  = model_vs[:-1][:][:]
    lvs  = model_vs[1:][:][:]
    urho = model_rho[:-1][:][:]
    lrho = model_rho[1:][:][:]
    rc_near=akirichards(uvp,uvs,urho,lvp,lvs,lrho,ang[0])
    rc_mid=akirichards(uvp,uvs,urho,lvp,lvs,lrho,ang[1])
    rc_far=akirichards(uvp,uvs,urho,lvp,lvs,lrho,ang[2])
    rc_near=np.concatenate((rc_near,np.zeros((1,n_traces))))  # add 1 row of zeros at the end
    rc_mid=np.concatenate((rc_mid,np.zeros((1,n_traces))))
    rc_far=np.concatenate((rc_far,np.zeros((1,n_traces))))
    return rc_near, rc_mid, rc_far


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def make_synth(rc,wavelet):
    '''
    Convolves reflectivities with wavelet.

    INPUT
    rc: 2D numpy array containing reflectivities
    wavelet

    OUTPUT
    synth: 2D numpy array containing seismic data

    Works with 1D arrays now (2015-05-07).
    '''
    nt=np.size(wavelet)
    if rc.ndim>1:
        [n_samples, n_traces] = rc.shape
        synth = np.zeros((n_samples+nt-1, n_traces))
        for i in range(n_traces):
            synth[:,i] = np.convolve(rc[:,i], wavelet)
        synth = synth[np.ceil(len(wavelet))/2:-np.ceil(len(wavelet))/2, :]
        synth=np.concatenate((synth,np.zeros((1,n_traces))))
    else:
        n_samples = rc.size
        synth = np.zeros(n_samples+nt-1)
        synth = np.convolve(rc, wavelet)
        synth = synth[np.ceil(len(wavelet))/2:-np.ceil(len(wavelet))/2]
        synth=np.concatenate((synth,[0]))
    return synth

# def make_synth(rc,wavelet):
#     nt=np.size(wavelet)
#     [n_samples, n_traces] = rc.shape
#     synth = np.zeros((n_samples+nt-1, n_traces))
#     for i in range(n_traces):
#         synth[:,i] = np.convolve(rc[:,i], wavelet)
#     synth = synth[np.ceil(len(wavelet))/2:-np.ceil(len(wavelet))/2, :]
#     synth=np.concatenate((synth,np.zeros((1,n_traces))))
#     return synth
#

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def make_synth_v2(rc,wavelet):
    '''
    Convolves reflectivities with wavelet.
    Alternative version using numpy apply_along_axis,
    slower than np.convolve with for loop.

    INPUT
    rc: 2D numpy array containing reflectivities
    wavelet

    OUTPUT
    synth: 2D numpy array containing seismic data
    '''
    nt=np.size(wavelet)
    [n_samples, n_traces] = rc.shape
    synth=np.zeros((n_samples+nt-1, n_traces))
    synth=np.apply_along_axis(lambda m: np.convolve(m,wavelet),axis=0,arr=rc)
    return synth

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def make_synth_v3(rc,wavelet):
    '''
    Convolves reflectivities with wavelet.
    Alternative version using scipy.ndimage.filters.convolve1d,
    slower than np.convolve with for loop.

    INPUT
    rc: 2D numpy array containing reflectivities
    wavelet

    OUTPUT
    synth: 2D numpy array containing seismic data
    '''
    from scipy.ndimage.filters import convolve1d
    nt=np.size(wavelet)
    [n_samples, n_traces] = rc.shape
    synth=np.zeros((n_samples+nt-1, n_traces))
    synth=convolve1d(rc,wavelet,axis=0)
    return synth

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def forward_model(model,aiprop,wavelet,dz,dt):
    """
    Meta function to do everything from scratch (zero-offset model).
    """
    earth = assign_ai(model, aiprop)
    vels = assign_vel(model, aiprop)
    earth_time=agilegeo.avo.depth_to_time(earth,vels,dz,dt,twt=True)
    rc = make_rc(earth_time)
    return make_synth(rc,wavelet)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def forward_model_elastic(model,elprop,wavelet,ang,dz,dt):
    """
    Meta function to do everything from scratch (angle-dependent models).
    """
    model_vp,model_vs,model_rho = assign_el(model,elprop)
    model_vp_time=agilegeo.avo.depth_to_time(model_vp,model_vp,dz,dt,twt=True)
    model_vs_time=agilegeo.avo.depth_to_time(model_vs,model_vp,dz,dt,twt=True)
    model_rho_time=agilegeo.avo.depth_to_time(model_rho,model_vp,dz,dt,twt=True)

    rc_near, rc_mid, rc_far=make_rc_elastic(model_vp_time,model_vs_time,model_rho_time,ang)
    near = make_synth(rc_near,wavelet)
    mid = make_synth(rc_mid,wavelet)
    far = make_synth(rc_far,wavelet)
    return near,mid,far

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def forward_model_elastic_decay(model,elprop,wav_near,wav_mid,wav_far,dz,dt):
    """
    Meta function to do everything from scratch (angle-dependent models).
    Uses angle-dependent wavelet to simulate frequency decay with offset.
    """
    model_vp,model_vs,model_rho = assign_el(model,elprop)
    model_vp_time=agilegeo.avo.depth_to_time(model_vp,model_vp,dz,dt,twt=True)
    model_vs_time=agilegeo.avo.depth_to_time(model_vs,model_vp,dz,dt,twt=True)
    model_rho_time=agilegeo.avo.depth_to_time(model_rho,model_vp,dz,dt,twt=True)

    rc_near, rc_mid, rc_far=make_rc_elastic(model_vp_time,model_vs_time,model_rho_time,ang)
    near = make_synth(rc_near,wav_near)
    mid = make_synth(rc_mid,wav_mid)
    far = make_synth(rc_far,wav_far)
    return near,mid,far

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def extract_amp(data,elprop,encasing_thickness,min_thickness,max_thickness,dt,freq):
    '''
    Extracts top and bottom real/apparent amplitudes from wedge.

    INPUT
    data: synthetic wedge in twt
    elprop: np.array([[vp1,rho1,vs1],[vp2,rho2,vs2],[vp3,rho3,vs3]])
    encasing_thickness
    min_thickness
    max_thickness
    dt: twt vertical sample rate

    OUTPUT
    toptwt0,bottwt0: top, bottom horizon (REAL)
    topamp0,botamp0: top, bottom amplitude (REAL)
    toptwt1,bottwt1: top, bottom horizon (APPARENT)
    topamp1,botamp1: top, bottom amplitude (APPARENT)
    '''
    [ns,nt]=data.shape
    twt=np.arange(0,ns*dt,dt)

    Fd=freq*1.3
    b=1/Fd
    cerca=int((b/dt)/2)

    # if Ip_above<Ip_below then we have an INCREASE in Ip = positive RC = peak
    top_is_peak=elprop[0,0]*elprop[0,1]<elprop[1,0]*elprop[1,1]
    bot_is_peak=elprop[1,0]*elprop[1,1]<elprop[2,0]*elprop[2,1]

    layer_1_twt=float(encasing_thickness)/elprop[0,0]*2
    incr=(max_thickness-min_thickness)/float(nt)

    toptwt0=np.zeros(nt)+layer_1_twt
    bottwt0=np.zeros(nt)+layer_1_twt+(min_thickness/elprop[1,0]*2)
    for i in range(nt):
        bottwt0[i]=bottwt0[i]+incr*i/elprop[1,0]*2

    # amplitude extraction at top,bottom REAL
    topamp0=np.zeros(nt)
    botamp0=np.zeros(nt)

    for i,val in enumerate(toptwt0):
        dd=np.abs(twt-val).argmin()
        window=data[dd,i]
        if top_is_peak:
            topamp0[i]=window.max()
        else:
            topamp0[i]=window.min()

    for i,val in enumerate(bottwt0):
        dd=np.abs(twt-val).argmin()
        window=data[dd,i]
        if bot_is_peak:
            botamp0[i]=window.max()
        else:
            botamp0[i]=window.min()

    # amplitude extraction at top,bottom APPARENT
    toptwt1=np.copy(toptwt0)
    bottwt1=np.copy(bottwt0)
    topamp1=np.zeros(nt)
    botamp1=np.zeros(nt)

    for i,val in enumerate(toptwt0):
        dd=np.abs(twt-val).argmin() # sample corresponding to horizon pick
        window=data[dd-cerca:dd+cerca,i] # amplitudes within a window centered on horizon pick and spanning -/+ samples (`cerca`)
        if np.any(window):
            if top_is_peak:
                toptwt1[i]=twt[np.abs(data[:,i]-window.max()).argmin()]
                topamp1[i]=window.max()
            else:
                toptwt1[i]=twt[np.abs(data[:,i]-window.min()).argmin()]
                topamp1[i]=window.min()
        else:
            toptwt1[i]=np.NaN
            topamp1[i]=np.NaN

    for i,val in enumerate(bottwt0):
        dd=np.abs(twt-val).argmin()
        window=data[dd-cerca:dd+cerca,i]
        if np.any(window):
            if bot_is_peak:
                bottwt1[i]=twt[np.abs(data[:,i]-window.max()).argmin()]
                botamp1[i]=window.max()
            else:
                bottwt1[i]=twt[np.abs(data[:,i]-window.min()).argmin()]
                botamp1[i]=window.min()
        else:
            bottwt1[i]=np.NaN
            botamp1[i]=np.NaN

    return toptwt0,bottwt0,topamp0,botamp0,toptwt1,bottwt1,topamp1,botamp1


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def extract_peakfreqs(data,min_thickness,max_thickness,dt):
    '''
    Extracts peak frequencies from wedge.

    INPUT
    data: synthetic wedge in twt
    min_thickness
    max_thickness
    dt: twt vertical sample rate

    OUTPUT
    aft: array with peak amplitude (A) at row 0, peak frequency (F) at row 1, thickness (T) at row 2
    spectra: array with amplitude spectra for all traces
    '''

    import aaplot
    from scipy.signal import argrelmax
    [ns,nt]=data.shape

    amp0,ff0=aaplot.ampspec(data[:,0],dt)
    spectra=np.zeros((amp0.size,nt))
    aft=np.zeros((3,nt)) # row 0: peak Amplitudes, row 1: peak Frequencies, row 2: Thickness
    for i in range(nt):
        amp,ff=aaplot.ampspec(data[:,i],dt)
        spectra[:,i]=amp
        peak_freq_list=ff[argrelmax(amp)]
        peak_amp_list=amp[argrelmax(amp)]
        if peak_freq_list.size==0:
            aft[0,i]=np.NaN
            aft[1,i]=np.NaN
        else:
            uu=peak_amp_list==np.max(peak_amp_list)
            peak_amp=peak_amp_list[uu]
            peak_freq=peak_freq_list[uu]
            aft[0,i]=peak_amp
            aft[1,i]=peak_freq
        incr=(max_thickness-min_thickness)/float(nt)
        aft[2,i]=i*incr+min_thickness
        # print peak_freq_list, peak_amp_list
        # print 'traccia %d, peak freq=%.2f, spessore=%.2f' % (i, peak_freq, ss[2,i])
    return aft, spectra
