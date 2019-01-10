import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# original data from Avseth, P., Mukerji, T. & Mavko, G. Quantitative Seismic Interpretation. (Cambridge University Press, 2005).
# available here:
# https://pangea.stanford.edu/researchgroups/srb/resources/books/quantitative-seismic-interpretation

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
import matplotlib.colors as colors
#            grey70     blue      red         brown
colori4 = ['#B3B3B3', '#003EFF','#FF0000',  '#996633']
cmap4   = colors.ListedColormap(colori4[0:len(colori4)], 'indexed')

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def colorbar_index(ncolors, cmap):
    mappable = plt.cm.ScalarMappable(cmap=cmap)
    mappable.set_array([])
    mappable.set_clim(-0.5, ncolors+0.5)
    colorbar = plt.colorbar(mappable)
    colorbar.set_ticks(np.linspace(0, ncolors, ncolors))
    colorbar.set_ticklabels(range(0,ncolors+1))
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


# Will use QSI dataset
# data for 5 wells are available, I will use well 2

w2=pd.read_table('well_2.txt', sep='\s+', header=None, skiprows=1, names=['DEPTH','VP','VS','RHO','GR', 'NPHI'], na_values=['-999.2500'])
tmp1=pd.read_table('well_2_denscorr.txt', sep='\s+', header=None, skiprows=1, names=['DEPTH','RHO_CORR'], na_values=['-999.2500'])
tmp2=pd.read_table('well_2_sats.txt', sep='\s+', header=None, skiprows=1, names=['DEPTH','SW','SWX'], na_values=['-999.2500'])
tmp2.DEPTH+=25  # adds KB?, i.e. after comparison with SW plot on p.262 (QSI) I see ~25m difference in SW log plot
                # so I reckon the depths in well_2_sats.txt must be TVDSS

#------------------------------------------------
# plot qc n.1: density vs corrected density
plt.figure()
w2.plot('RHO', 'DEPTH', style=':r',xlim=(1.5, 3.0))  # will be renamed RHO_OLD
tmp1.plot('RHO_CORR', 'DEPTH', style='-k',xlim=(1.5, 3.0))
plt.gca().invert_yaxis()

plt.figure()
w2.plot('VP', 'DEPTH', style='-r')
plt.xlim(1000,4000), plt.ylim(2140,2180)
plt.gca().invert_yaxis()


# plot qc n.2: shifted SW log
plt.figure(); tmp2.plot('SW', 'DEPTH', style='-b',xlim=(-.1, 1.1), ylim=(2000, 2600)); plt.gca().invert_yaxis()
#------------------------------------------------

# will import RHO_CORR as the good one, therefore current RHO will be renamed to RHO_OLD
w2.rename(columns={'RHO':'RHO_OLD'}, inplace=True)

# add RHO_CORR as RHO to DataFrame w2 using np.interp
w2['RHO']= np.interp(w2.DEPTH.values, tmp1.DEPTH.values, tmp1.RHO_CORR.values, left=np.NaN, right=np.NaN)

# add SW, SWX to DataFrame w2 using np.interp
w2['SW']=  np.interp(w2.DEPTH.values, tmp2.DEPTH.values, tmp2.SW.values,  left=np.NaN, right=np.NaN)
w2['SWX']= np.interp(w2.DEPTH.values, tmp2.DEPTH.values, tmp2.SWX.values, left=np.NaN, right=np.NaN)

# calculates Ip, Is, Vp/Vs and transform vels in m/s
w2.VP=w2.VP*1000
w2.VS=w2.VS*1000
w2['VPVS']=w2.VP/w2.VS
w2['IP']=w2.VP*w2.RHO
w2['IS']=w2.VS*w2.RHO

# input elastic parameters; see also QSI, p.261, 336, 338
rho_qz=2.65;  k_qz=37;  mu_qz=44
rho_sh=2.81;  k_sh=15;  mu_sh=5
rho_w=1.09;   k_w=2.8
rho_o=0.78;   k_o=0.94 # # oil gravity: 32 API, GOR:    64


# shale volume (%) = (GR - GRmin)/(GRmax - GRmin)  QSI reports WRONG equation!!! (GRmax-GR)/(GRmax-Grmin)
# see here http://www.spec2000.net/11-vshgr.htm
w2['VSH']=(w2.GR-w2.GR.min())/(w2.GR.max()-w2.GR.min())

# rho_matrix = vol_min1*rho_min1 + vol_min2*rho_min2
w2['RHOm']=w2.VSH*rho_sh + (1-w2.VSH)*rho_qz

# rho_fluid = Sw*rho_water + (1-Sw)*rho_oil
w2['RHOf']=w2.SW*rho_w + (1-w2.SW)*rho_o

# porosity=(rho_matrix- rho_log)/(rho_matrix - rho_fluid)
# rho                   = (1-phi)*rho_m     + phi*rho_f
# rho                   = rho_m - phi*rho_m + phi*rho_f
# phi*rho_m - phi*rho_f = (rho_m - rho)
# phi*(rho_m - rho_f)   = (rho_m - rho)
# phi                   = (rho_m - rho) / (rho_m - rho_f)
w2['PHI']= (w2.RHOm-w2.RHO) / ( w2.RHOm- w2.RHOf)


w2.to_csv('qsiwell2.csv',index=False)

# ********************************************************************
# ********************************************************************
# ********************************************************************

w1=pd.read_table('well_1.txt', sep='\s+', header=None, skiprows=1, names=['DEPTH','VP','RHO','GR'])
w5=pd.read_table('well_5.txt', sep='\s+', header=None, skiprows=1, names=['DEPTH','DTP','DTS','GR', 'RHO'])

# ********************************************************************
# ********************************************************************
# ********************************************************************

# log plot simile a quello che si ottiene con aageofisica.wellplot, con in piu' log di facies
ll=w2
versionelitolog=1
ztop=2100
zbot=2400
dummy=np.zeros(len(ll.VCL))


velmin=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['VS']].min().values
velmax=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['VP']].max().values
ipmin=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['IP']].min().values
ipmax=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['IP']].max().values
rmin=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['VPVS']].min().values
rmax=ll[(ll.DEPTH>=ztop) & (ll.DEPTH<=zbot)].ix[:,['VPVS']].max().values

f, ax = plt.subplots(nrows=1, ncols=5, sharey=True, figsize=(12,6))
ll.plot(x='SW',   y='DEPTH', ax=ax[0], style='b', label='Sw');
ll.plot(x='VCL',  y='DEPTH', ax=ax[0], style='g', label='Vcl');
ll.plot(x='PHI',  y='DEPTH', ax=ax[0], style='k', label='phi');
ll.plot(x='VP',   y='DEPTH', ax=ax[1], style='k');
ll.plot(x='VS',   y='DEPTH', ax=ax[1], style='r');
ll.plot(x='IP',   y='DEPTH', ax=ax[2], style='k');
ll.plot(x='VPVS', y='DEPTH', ax=ax[3], style='k');
if versionelitolog==0:
    ax[4].plot(dummy[(ll.LFC==0).values],ll.DEPTH[ll.LFC==0],'s',color='#B3B3B3',markeredgewidth=0)
    ax[4].plot(dummy[(ll.LFC==1).values],ll.DEPTH[ll.LFC==1],'sb',label='LFC 1',markeredgewidth=0)
    ax[4].plot(dummy[(ll.LFC==2).values],ll.DEPTH[ll.LFC==2],'sr',label='LFC 2',markeredgewidth=0)
    ax[4].plot(dummy[(ll.LFC==3).values],ll.DEPTH[ll.LFC==3],'s',color='#996633',label='LFC 3',markeredgewidth=0)
    ax[4].set_xlabel('LFC'),              ax[4].set_xlim(-0.5,3),       ax[3].set_ylim(ztop,zbot)
else:
    ax[4].plot(ll.LFC[ll.LFC==0],ll.DEPTH[ll.LFC==0],'s',color='#B3B3B3',markeredgewidth=0)
    ax[4].plot(ll.LFC[ll.LFC==1],ll.DEPTH[ll.LFC==1],'sb',label='LFC 1',markeredgewidth=0)
    ax[4].plot(ll.LFC[ll.LFC==2],ll.DEPTH[ll.LFC==2],'sr',label='LFC 2',markeredgewidth=0)
    ax[4].plot(ll.LFC[ll.LFC==3],ll.DEPTH[ll.LFC==3],'s',color='#996633',label='LFC 3',markeredgewidth=0)
    ax[4].set_xlabel('LFC'),              ax[4].set_xlim(-1,10),       ax[3].set_ylim(ztop,zbot)
ax[0].set_xlabel('Vcl/Sw/phi'),       ax[0].set_xlim(-0.1,1.1),                          ax[0].set_ylim(ztop,zbot)
ax[1].set_xlabel('Velocities (m/s)'), ax[1].set_xlim(velmin-velmin*.1,velmax+velmax*.1), ax[1].set_ylim(ztop,zbot)
ax[2].set_xlabel('Ip (m/s*g/cc)'),    ax[2].set_xlim(ipmin-ipmin*.1,ipmax+ipmax*.1),     ax[2].set_ylim(ztop,zbot)
ax[3].set_xlabel('Vp/Vs'),            ax[3].set_xlim(rmin-rmin*.01,rmax+rmax*.01),       ax[3].set_ylim(ztop,zbot)
ax[0].invert_yaxis()
ax[0].legend(fontsize='small', loc='lower left')
ax[1].legend(fontsize='small', loc='lower left')
ax[4].legend(fontsize='small', loc='lower right')
ax[4].set_xticklabels([])