import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def gassmann(vp1, vs1, rho1, rhof1, kfluid1, rhof2, kfluid2, kmin, phi):
    '''
    Calculate elastic properties (Vp, Vs, density) after
    replacing fluids via Gassmann's equation.
    (aadm 2015-2018)

    Parameters
    ----------
    vp1, vs1 : float or array_like
        Initial P-wave and S-wave velocity in m/s.
    rho1 : float or array_like
        Final density in g/cm3.
    rhof1 : float or array_like
        Initial fluid density in g/cm3.
    kfluid1 : float or array_like
        Initial fluid bulk modulus in GPa.
    rhof2 : float or array_like
        Final fluid density in g/cm3.
    kfluid2 : float or array_like
        Final fluid bulk modulus in GPa.
    kmin : float or array_like
        Mineral bulk modulus in GPa.
    phi : float or array_like
        Porosity in fraction.

    Returns
    -------
    vp2, vs2 : float or array_like
        Final P-wave and S-wave velocity in m/s.
    rho2 : float or array_like
        Final density in g/cm3.
    ksat2 : float or array_like
        Final saturated rock bulk modulus in GPa.
    kdry : float or array_like
        Dry-rock bulk modulus in GPa.
    '''
    # convert density to kg/m3 and elastic moduli to Pa
    d1 = rho1*1e3
    df1 = rhof1*1e3
    df2 = rhof2*1e3
    k0 = kmin*1e9
    kf1 = kfluid1*1e9
    kf2 = kfluid2*1e9
    d2 = d1 - phi * df1 + phi * df2
    mu1 = d1 * vs1**2
    k1 = d1 * vp1**2 - (4/3) * mu1
    kd = (k1 * (phi*k0 / kf1 + 1-phi) - k0) / (phi*k0 / kf1 + k1 / k0 - 1-phi)
    mu2 = mu1
    with np.errstate(divide='ignore', invalid='ignore'):
        k2 = kd + (1 - kd/k0 )**2 / (phi/kf2 + (1-phi)/k0 - kd/k0**2)
        vp2 = np.sqrt((k2 + 4/3*mu2) / d2)
        vs2 = np.sqrt(mu2 / d2)
        return vp2, vs2, d2/1e3, k2/1e9, kd/1e9

    
def bulk(vp, vs, rho):
    '''
    Calculate bulk modulus.
    (aadm 2017)

    Parameters
    ----------
    vp, vs: float
        P- and S-wave velocity [m/s].
    rho: float
        Density [g/cm3].

    Returns
    -------
    K: float
        Bulk modulus [GPa].
    '''
    # converts density to SI (kg/m3)
    D = rho*1e3
    K = D*vp**2 - 4/3*D*vs**2
    return K/1e9


def shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta, approx=True, terms=False):
    '''
    Calculate P-wave reflectivity with Shuey's equation.
    (aadm 2016, 2020)

    Parameters
    ----------
    vp1, vs1, rho1 : float or array_like
        P-, S-wave velocity (m/s) and density (g/cm3) of upper medium.
    vp2, vs2, rho2 : float or array_like
        P-, S-wave velocity (m/s) and density (g/cm3) of lower medium.
    theta : int or array_like
        Angle of incidence (degree).
    approx : bool, optional
        If True returns approximate 2-terms form. Default: True
    terms :  bool, optional
        If True returns reflectivity, intercept and gradient.
        Default: False.

    Returns
    -------
    R : float or array_like
        Reflectivity at angle theta.
    R0, G : float
        Intercept and gradient, only output if terms is True.

    Notes
    -----
    If input properties are arrays with length n and angles are also arrays with length m, the function returns a (n,m) array.
    
    References
    ----------
    Avseth et al. (2005), Quantitative Seismic Interpretation, Cambridge University Press (p.182)
    '''
    a = np.radians(theta)
    dvp = vp2-vp1
    dvs = vs2-vs1
    drho = rho2-rho1
    vp  = np.mean([vp1, vp2], axis=0)
    vs  = np.mean([vs1, vs2], axis=0)
    rho = np.mean([rho1, rho2], axis=0)
    R0 = 0.5*(dvp/vp + drho/rho)
    G  = 0.5*(dvp/vp) - 2*(vs**2/vp**2)*(drho/rho+2*(dvs/vs))
    F =  0.5*(dvp/vp)
    # if angles is an array
    if a.size>1:
        R0 = R0.reshape(-1,1)
        G = G.reshape(-1,1)
        F = F.reshape(-1,1)
    if approx:
        R = R0 + G*np.sin(a)**2
    else:
        R = R0 + G*np.sin(a)**2 + F*(np.tan(a)**2-np.sin(a)**2)
    if terms:
        return R, R0, G
    else:
        return R
    
    
def gassmann_phi(ksat1, ksat2, kf1, kf2, kmin):
    a = (kmin-ksat1)*(kmin-ksat2)*(kf1-kf2) 
    b = (kmin-kf1)*(kmin-kf2)*(ksat1-ksat2)
    return a / b


def percdiff(start, end):
    return (end-start)/start


def classref(near=5, far=30, above=None, below=None, mx=.6, plot_brine=True):
    tmp_shl = np.array([[3094, 1515, 2.40, 0],
                        [2643, 1167, 2.29, 0],
                        [2192, 818, 2.16, 0],
                        [3240, 1620, 2.34, 0]])
    tmp_ssg = np.array([[4050, 2526, 2.21, .2],
                        [2781, 1665, 2.08, .25],
                        [1542, 901, 1.88, .33],
                        [1650, 1090, 2.07, .156]])
    tmp_ssb = np.array([[4115, 2453, 2.32, .2],
                        [3048, 1595, 2.23, .25],
                        [2134, 860, 2.11, .33],
                        [2590, 1060, 2.21, .156]])
    avocl = ['CLASS-1', 'CLASS-2', 'CLASS-3', 'CLASS-4']
    logs = ['VP', 'VS', 'RHO', 'PHI']
    shl = pd.DataFrame(tmp_shl, columns=logs, index=avocl)
    ssg = pd.DataFrame(tmp_ssg, columns=logs, index=avocl)
    ssb = pd.DataFrame(tmp_ssb, columns=logs, index=avocl)

    opttxt = dict(weight='bold', ha='left', va='center')
    mrkg = {'ms': 10, 'mew': 2, 'ls': 'none'}
    mrkb = {'ms': 10, 'mew': 2, 'ls': 'none', 'mfc':'none'}
    mrk_sel = {'marker': '*', 'mec': 'k', 'mfc': 'white', 'ms': 16, 'ls': 'none', 'mew': 2}

    angs = np.array([near, far])
    tmp = ['C0', 'C1', 'C2', 'C3']
    cc = dict(zip(avocl, tmp))
    tmp = ['s', 'P', 'v', '^']
    mm  = dict(zip(avocl, tmp))

    f, ax = plt.subplots(constrained_layout=True)
    ax.axhline(0, color='k', lw=3)
    ax.axvline(0, color='k', lw=3)
    for i, sh in shl.iterrows():
        vpsh, vssh, dsh = sh['VP'], sh['VS'], sh['RHO']
        vpb, vsb, db = ssb.loc[i, 'VP'], ssb.loc[i, 'VS'], ssb.loc[i, 'RHO']
        vpg, vsg, dg = ssg.loc[i, 'VP'], ssg.loc[i, 'VS'], ssg.loc[i, 'RHO']
        Ab, Ib, Gb = shuey(vpsh, vssh, dsh, vpb, vsb, db, angs, terms=True)
        Ag, Ig, Gg = shuey(vpsh, vssh, dsh, vpg, vsg, dg, angs, terms=True)
        ax.plot(Ig, Gg, fillstyle='full', label=sh.name, marker=mm[i], mfc=cc[i], mec=cc[i], **mrkg)
        if plot_brine:
            ax.plot(Ib, Gb, fillstyle='none', label=None, marker=mm[i], mec=cc[i], **mrkb)
    if above is not None:
        vp0, vs0, d0 = above
        vp1, vs1, d1 = below
        Ax, Ix, Gx = shuey(vp0, vs0, d0, vp1, vs1, d1, angs, terms=True)
        ax.plot(Ix, Gx, **mrk_sel)
    ax.set_xlabel('Intercept')
    ax.set_ylabel('Gradient')
    ax.legend()
    ax.set_xlim(-mx, mx)
    ax.set_ylim(-mx, mx)
    ax.set_aspect('equal', 'box')
    ax.grid()
    if plot_brine:
        ax.set_title('Filled markers: gas, empty=brine')