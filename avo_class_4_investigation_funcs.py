import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def percdiff(start, end):
    return (end-start)/start


def find_nearest(a, a0):
    '''
    Element in nd array `a` closest to the scalar value `a0`
    '''
    idx = np.abs(a - a0).argmin()
    return idx, a[idx]


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


def gassmann_phi(ksat1, ksat2, kf1, kf2, kmin):
    a = (kmin-ksat1)*(kmin-ksat2)*(kf1-kf2) 
    b = (kmin-kf1)*(kmin-kf2)*(ksat1-ksat2)
    return a / b


def hertzmindlin(K0, G0, sigma, phi_c=0.4, Cn=8.6, f=1):
    '''
    Hertz-Mindlin model.
    (aadm 2015)

    Parameters
    ----------   
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in GPa.
    phi : float or array_like
        Porosity.
    sigma : float
        Effective stress in MPa.
    phi_c : float, optional
        Critical porosity. Default: 0.4
    Cn : float, optional
        Coordination number Default: 8.6.
    f : float, optional
        Shear modulus correction factor,
        f=1 for dry pack with perfect adhesion
        between particles and f=0 for dry frictionless pack.

    Returns
    -------
    Kdry, Gdry : float or array_like
        Dry rock bulk & shear modulus in GPa.

    References
    ----------
    Mavko et al. (2009), The Rock Physics Handbook, Cambridge University Press (p.246)
    '''
    sigma0 = sigma / 1e3  # converts pressure in same units as solid moduli (GPa)
    pr0 = (3*K0-2*G0) / (6*K0+2*G0)  # poisson's ratio of mineral mixture
    Khm = (sigma0*(Cn**2*(1 - phi_c)**2*G0**2) / (
           18*np.pi**2 * (1 - pr0)**2))**(1/3)
    Ghm = ((2+3*f-pr0*(1+3*f)) / (5*(2-pr0))) * (
          (sigma0 * (3 * Cn**2 * (1 - phi_c)**2 * G0**2) / (
           2 * np.pi**2 * (1 - pr0)**2)))**(1/3)
    return Khm, Ghm


def softsand(K0, G0, phi, sigma, phi_c=0.4, Cn=8.6, f=1):
    '''
    Soft sand, or friable sand or uncemented sand model.
    (aadm 2015)

    Parameters
    ----------
    K0, G0 : float or array_like
        Mineral bulk & shear modulus in GPa.
    phi : float or array_like
        Porosity.
    sigma : float
        Effective stress in MPa.
    phi_c : float, optional
        Critical porosity. Default: 0.4
    Cn : float, optional
        Coordination number Default: 8.6.
    f : float, optional
        Shear modulus correction factor,
        f=1 for dry pack with perfect adhesion
        between particles and f=0 for dry frictionless pack.

    Returns
    -------
    Kdry, Gdry : float or array_like
        Dry rock bulk & shear modulus in GPa.

    References
    ----------
    Mavko et al. (2009), The Rock Physics Handbook, Cambridge University Press (p.258)
    '''
    Khm, Ghm = hertzmindlin(K0, G0, sigma, phi_c, Cn, f)
    Kdry = -4/3 * Ghm + (((phi / phi_c) / (Khm + 4/3 * Ghm)) + (
           (1 - phi / phi_c) / (K0 + 4/3 * Ghm)))**-1
    gxx = Ghm / 6 * ((9 * Khm + 8 * Ghm) / (Khm + 2 * Ghm))
    Gdry = -gxx + ((phi / phi_c) / (Ghm + gxx) + (
           (1 - phi / phi_c) / (G0 + gxx)))**-1
    return Kdry, Gdry


def vels(kdry, gdry, kmin, rho0, kfluid, rhof, phi):
    '''
    Calculate velocities and densities of saturated rock
    using Gassmann's equation.
    (aadm 2015-2018)

    Parameters
    ----------
    kdry : float or array_like
        Dry-rock bulk modulus in GPa.
    gdry : float or array_like
        Dry-rock shear modulus in GPa.
    kmin : float or array_like
        Mineral bulk modulus in GPa.
    rho0 : float or array_like
        Mineral density in g/cm3.
    kfluid : float or array_like
        Fluid bulk modulus in GPa.
    rhof : float or array_like
        Fluid density in g/cm3.
    phi : float or array_like
        Porosity in fraction.

    Returns
    -------
    vp, vs : float or array_like
        P-wave and S-wave velocity in m/s.
    rho : float or array_like
        Density in g/cm3.
    ksat : float or array_like
        Saturated rock bulk modulus in GPa.
    '''
    # convert density to kg/m3 and elastic moduli to Pa
    d0 = rho0*1e3
    df = rhof*1e3
    kd = kdry*1e9
    gd = gdry*1e9
    kf = kfluid*1e9
    k0 = kmin*1e9
    rho = d0 * (1 - phi) + df * phi
    with np.errstate(divide='ignore', invalid='ignore'):
        ksat = kd + (1 - kd/k0)**2 / (phi/kf + (1-phi)/k0 - kd/k0**2)
        vp = np.sqrt((ksat+4/3*gd)/rho)
        vs = np.sqrt(gd/rho)
    return vp, vs, rho/1e3, ksat/1e9


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


def classref(near=5, far=30, mx=.6, plot_brine=False, plot_colorzones=True):
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
    avocl = ['CLASS1', 'CLASS2', 'CLASS3', 'CLASS4']
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
    ax.set_xlabel('Intercept')
    ax.set_ylabel('Gradient')
    ax.legend()
    ax.set_xlim(-mx, mx)
    ax.set_ylim(-mx, mx)
    ax.set_aspect('equal', 'box')
    ax.grid()
    if plot_colorzones:
        opt1 = dict(edgecolor='None', alpha=0.2)
        cl1_area = patches.Rectangle((0.02, -1), .98, 1, facecolor=cc['CLASS1'], **opt1)
        cl2_area = patches.Rectangle((-0.02, -1), .04, 2, facecolor=cc['CLASS2'], **opt1)
        cl3_area = patches.Rectangle((-1, -1), .98, 1, facecolor=cc['CLASS3'], **opt1)
        cl4_area = patches.Rectangle((-1, 0), .98, 1, facecolor=cc['CLASS4'], **opt1)
        background = patches.Polygon([[-1, 1], [1, -1], [1, 1]], facecolor='w')
        ax.add_patch(cl1_area)
        ax.add_patch(cl2_area)
        ax.add_patch(cl3_area)
        ax.add_patch(cl4_area)
        ax.add_patch(background)
        opt2 = dict(ha='center', va='center', weight='bold', size='large')
        ax.text(.15, -.3, 'Class 1', color=cc['CLASS1'], **opt2)
        ax.text(0, -.25, 'Class 2/2p', color=cc['CLASS2'], **opt2)
        ax.text(-.35, -.3, 'Class 3', color=cc['CLASS3'], **opt2)
        ax.text(-.35, .15, 'Class 4', color=cc['CLASS4'], **opt2)
    if plot_brine:
        ax.set_title('Filled markers: gas, empty=brine')
    return ax