# from https://sites.google.com/site/theodoregoetz/notes/matplotlib_colormapadjust

import math, copy
import numpy
import numpy as np
from matplotlib import pyplot, colors, cm
import matplotlib
from scipy import interpolate

def cmap_discretize(cmap, N):
    """Return a discrete colormap from the continuous colormap cmap.

        cmap: colormap instance, eg. cm.jet.
        N: Number of colors.

    Example
        x = resize(arange(100), (5,100))
        djet = cmap_discretize(cm.jet, 5)
        imshow(x, cmap=djet)
    """

    cdict = cmap._segmentdata.copy()
    # N colors
    colors_i = np.linspace(0,1.,N)
    # N+1 indices
    indices = np.linspace(0,1.,N+1)
    for key in ('red','green','blue'):
        # Find the N colors
        D = np.array(cdict[key])
        I = interpolate.interp1d(D[:,0], D[:,1])
        colors = I(colors_i)
        # Place these colors at the correct indices.
        A = np.zeros((N+1,3), float)
        A[:,0] = indices
        A[1:,1] = colors
        A[:-1,2] = colors
        # Create a tuple for the dictionary.
        L = []
        for l in A:
            L.append(tuple(l))
        cdict[key] = tuple(L)
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)


def cmap_powerlaw_adjust(cmap, a):
    '''
    returns a new colormap based on the one given
    but adjusted via power-law:

    newcmap = oldcmap**a
    '''
    if a < 0.:
        return cmap
    cdict = copy.copy(cmap._segmentdata)
    fn = lambda x : (x[0]**a, x[1], x[2])
    for key in ('red','green','blue'):
        cdict[key] = map(fn, cdict[key])
        cdict[key].sort()
        assert (cdict[key][0]<0 or cdict[key][-1]>1), \
            "Resulting indices extend out of the [0, 1] segment."
    return colors.LinearSegmentedColormap('colormap',cdict,1024)

def cmap_center_adjust(cmap, center_ratio):
    '''
    returns a new colormap based on the one given
    but adjusted so that the old center point higher
    (>0.5) or lower (<0.5)
    '''
    if not (0. < center_ratio) & (center_ratio < 1.):
        return cmap
    a = math.log(center_ratio) / math.log(0.5)
    return cmap_powerlaw_adjust(cmap, a)

def cmap_center_point_adjust(cmap, range, center):
    '''
    converts center to a ratio between 0 and 1 of the
    range given and calls cmap_center_adjust(). returns
    a new adjusted colormap accordingly
    '''
    if not ((range[0] < center) and (center < range[1])):
        return cmap
    return cmap_center_adjust(cmap,
        abs(center - range[0]) / abs(range[1] - range[0]))


if __name__ == '__main__':
    ### create some 2D histogram-type data
    def func3(x,y):
        return (1- x/2 + x**5 + y**3)*numpy.exp(-x**2-y**2)
    x = numpy.linspace(-3.0, 3.0, 60)
    y = numpy.linspace(-3.0, 3.0, 60)
    X,Y = numpy.meshgrid(x, y)
    Z = func3(X, Y)
    extent = [x[0],x[-1],y[0],y[-1]]


    plotkwargs = {
        'extent' : extent,
        'origin' : 'lower',
        'interpolation' : 'nearest',
        'aspect' : 'auto'}

    fig = pyplot.figure(figsize=(8,3))
    fig.subplots_adjust(left=.05,bottom=.11,right=.94,top=.83,wspace=.35)
    ax = [fig.add_subplot(1,3,i) for i in range(1,4,1)]

    cmap = cm.seismic

    plt = ax[0].imshow(Z, cmap=cmap, **plotkwargs)
    cb = ax[0].figure.colorbar(plt, ax=ax[0])
    ax[0].set_title('cmap: seismic')

    plt = ax[1].imshow(Z, cmap=cmap_center_adjust(cmap, 0.75), **plotkwargs)
    cb = ax[1].figure.colorbar(plt, ax=ax[1])
    ax[1].set_title('center raised by 25%')

    plt = ax[2].imshow(Z,
        cmap=cmap_center_point_adjust(cmap,[numpy.min(Z),numpy.max(Z)],0),
        **plotkwargs)
    cb = ax[2].figure.colorbar(plt, ax=ax[2])
    ax[2].set_title('center set to 0')

    pyplot.show()
