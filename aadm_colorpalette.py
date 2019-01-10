import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import  LinearSegmentedColormap
import matplotlib.cm as cm
import os

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# https://mycarta.wordpress.com/color-palettes/
def creacolormap(arr):
    b3=arr[:,2] # value of blue at sample n
    b2=arr[:,2] # value of blue at sample n
    b1=np.linspace(0,1,len(b2)) # position of sample n - ranges from 0 to 1
    g3=arr[:,1]
    g2=arr[:,1]
    g1=np.linspace(0,1,len(g2))
    r3=arr[:,0]
    r2=arr[:,0]
    r1=np.linspace(0,1,len(r2))
    # creating list
    R=zip(r1,r2,r3)
    G=zip(g1,g2,g3)
    B=zip(b1,b2,b3)
    # transposing list
    RGB=zip(R,G,B)
    rgb=zip(*RGB)
    # creating dictionary
    k=['red', 'green', 'blue']
    colori = dict(zip(k,rgb)) # makes a dictionary from 2 lists
    return colori

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# http://schubert.atmos.colostate.edu/~cslocum/custom_cmap.html
def make_cmap(colors, position=None, bit=False):
    '''
    make_cmap takes a list of tuples which contain RGB values. The RGB
    values may either be in 8-bit [0 to 255] (in which bit must be set to
    True when called) or arithmetic [0 to 1] (default). make_cmap returns
    a cmap with equally spaced colors.
    Arrange your tuples so that the first color is the lowest value for the
    colorbar and the last is the highest.
    position contains values from 0 to 1 to dictate the location of each color.
    '''
    bit_rgb = np.linspace(0,1,256)
    if position == None:
        position = np.linspace(0,1,len(colors))
    else:
        if len(position) != len(colors):
            sys.exit("position length must be the same as colors")
        elif position[0] != 0 or position[-1] != 1:
            sys.exit("position must start with 0 and end with 1")
    if bit:
        for i in range(len(colors)):
            colors[i] = (bit_rgb[colors[i][0]],
                         bit_rgb[colors[i][1]],
                         bit_rgb[colors[i][2]])
    cdict = {'red':[], 'green':[], 'blue':[]}
    for pos, color in zip(position, colors):
        cdict['red'].append((pos, color[0], color[0]))
        cdict['green'].append((pos, color[1], color[1]))
        cdict['blue'].append((pos, color[2], color[2]))

    cmap = LinearSegmentedColormap('my_colormap',cdict,256)
    return cmap
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# http://stackoverflow.com/questions/3279560/invert-colormap-in-matplotlib
def reverse_colourmap(cmap, name = 'my_cmap_r'):
    """
    In:
    cmap, name
    Out:
    my_cmap_r

    Explanation:
    t[0] goes from 0 to 1
    row i:   x  y0  y1 -> t[0] t[1] t[2]
                   /
                  /
    row i+1: x  y0  y1 -> t[n] t[1] t[2]

    so the inverse should do the same:
    row i+1: x  y1  y0 -> 1-t[0] t[2] t[1]
                   /
                  /
    row i:   x  y1  y0 -> 1-t[n] t[2] t[1]
    """
    reverse = []
    k = []

    for key in cmap._segmentdata:
        k.append(key)
        channel = cmap._segmentdata[key]
        data = []

        for t in channel:
            data.append((1-t[0],t[2],t[1]))
        reverse.append(sorted(data))

    LinearL = dict(zip(k,reverse))
    my_cmap_r = LinearSegmentedColormap(name, LinearL)
    return my_cmap_r

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#  adattato da http://matplotlib.org/examples/color/colormaps_reference.html
def plot_color_gradients(colormaps):
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    nrows = len(colormaps)
    fig, axes = plt.subplots(nrows=nrows)
    fig.subplots_adjust(top=0.95, bottom=0.01, left=0.2, right=0.99)
    for ax, name in zip(axes, colormaps):
        ax.imshow(gradient, aspect='auto', cmap=name)
        pos = list(ax.get_position().bounds)
        x_text = pos[0]
        y_text = pos[1] + pos[3]/2.
        fig.text(x_text, y_text, name, va='center', ha='right', fontsize=10)
    for ax in axes:
        ax.set_axis_off()


#==== colorbars di Matteo Niccoli
cc = np.loadtxt('colormap_sawtooth.csv', delimiter=',')
cmap_jetsaw = LinearSegmentedColormap('jetsaw',creacolormap(cc))
cm.register_cmap(name='jetsaw', cmap=cmap_jetsaw)

cmap_jetsaw_r = LinearSegmentedColormap('jetsaw_r',creacolormap(np.flipud(cc)))
cm.register_cmap(name='jetsaw_r', cmap=cmap_jetsaw_r)

#==== colorbars di Peter Kovesi http://peterkovesi.com/projects/colourmaps/
path = 'CETperceptual_csv_0_1'
for ff in os.listdir(path):
    pippo=ff.split(sep='_')[0:3]
    cname='_'.join(pippo)
#    print('{:>50s} --> {:<30s}'.format(ff, cname))
    cc = np.loadtxt(os.path.join(path,ff), delimiter=',')
    cm.register_cmap(name=cname, cmap=LinearSegmentedColormap(cname,creacolormap(cc)))
    cm.register_cmap(name=cname+'_r', cmap=LinearSegmentedColormap(cname,creacolormap(np.flipud(cc))))



#==== colorbar di Decision Space
 # black  purple       blue       cyan         green      red       yellow
colors =[(0,0,0), (255,0,255), (0,0,255), (0,255,255), (0,255,0), (255,0,0), (255,255,0), (255,255,255)]
cmap_landmark = make_cmap(colors, bit=True)
cm.register_cmap(name='landmark', cmap=cmap_landmark)

cmap_landmark_r = reverse_colourmap(cmap_landmark)
cm.register_cmap(name='landmark_r', cmap=cmap_landmark_r)


#==== pulizia
del cc,colors,cmap_jetsaw,cmap_jetsaw_r,cmap_landmark,cmap_landmark_r
del cm, LinearSegmentedColormap, creacolormap, make_cmap, reverse_colourmap
