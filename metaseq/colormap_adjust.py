"""
Module to handle custom colormaps.

`cmap_powerlaw_adjust`, `cmap_center_adjust`, and
`cmap_center_adjust` are from
https://sites.google.com/site/theodoregoetz/notes/matplotlib_colormapadjust
"""

import math
import copy
import numpy
import numpy as np
from matplotlib import pyplot, colors, cm
import matplotlib
import colorsys


def color_test(color):
    """
    Figure filled in with `color`; useful for troubleshooting or experimenting
    with colors
    """
    if isinstance(color, np.ndarray):
        color = color.ravel()
    fig = pyplot.figure(figsize=(2, 2), facecolor=color)


def smart_colormap(vmin, vmax, color_high='#b11902', hue_low=0.6):
    """
    Creates a "smart" colormap that is centered on zero, and accounts for
    asymmetrical vmin and vmax by matching saturation/value of high and low
    colors.

    It works by first creating a colormap from white to `color_high`.  Setting
    this color to the max(abs([vmin, vmax])), it then determines what the color
    of min(abs([vmin, vmax])) should be on that scale.  Then it shifts the
    color to the new hue `hue_low`, and finally creates a new colormap with the
    new hue-shifted as the low, `color_high` as the max, and centered on zero.

    :param color_high: a matplotlib color -- try "#b11902" for a nice red
    :param hue_low: float in [0, 1] -- try 0.6 for a nice blue
    :param vmin: lowest value in data you'll be plotting
    :param vmax: highest value in data you'll be plotting
    """
    # first go from white to color_high
    orig_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'test', ['#FFFFFF', color_high], N=2048)

    # For example, say vmin=-3 and vmax=9.  If vmin were positive, what would
    # its color be?
    vmin = float(vmin)
    vmax = float(vmax)
    mx = max([vmin, vmax])
    mn = min([vmin, vmax])
    frac = abs(mn / mx)
    rgb = orig_cmap(frac)[:-1]

    # Convert to HSV and shift the hue
    hsv = list(colorsys.rgb_to_hsv(*rgb))
    hsv[0] = hue_low
    new_rgb = colorsys.hsv_to_rgb(*hsv)
    new_hex = matplotlib.colors.rgb2hex(new_rgb)

    zeropoint = abs(-vmin / (vmax - vmin))

    # Create a new colormap using the new hue-shifted color as the low end
    new_cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
        'test', [(0, new_rgb), (zeropoint, '#FFFFFF'), (1, color_high)],
        N=2048)

    return new_cmap


def cmap_powerlaw_adjust(cmap, a):
    """
    Returns a new colormap based on the one given
    but adjusted via power-law, `newcmap = oldcmap**a`.

    :param cmap: colormap instance (e.g., cm.jet)
    :param a: power
    """
    if a < 0.:
        return cmap
    cdict = copy.copy(cmap._segmentdata)
    fn = lambda x: (x[0] ** a, x[1], x[2])
    for key in ('red', 'green', 'blue'):
        cdict[key] = map(fn, cdict[key])
        cdict[key].sort()
        assert (cdict[key][0] < 0 or cdict[key][-1] > 1), \
            "Resulting indices extend out of the [0, 1] segment."
    return colors.LinearSegmentedColormap('colormap', cdict, 1024)


def cmap_center_adjust(cmap, center_ratio):
    """
    Returns a new colormap based on the one given
    but adjusted so that the old center point higher
    (>0.5) or lower (<0.5)

    :param cmap: colormap instance (e.g., cm.jet)
    :param center_ratio:

    """
    if not (0. < center_ratio) & (center_ratio < 1.):
        return cmap
    a = math.log(center_ratio) / math.log(0.5)
    return cmap_powerlaw_adjust(cmap, a)


def cmap_center_point_adjust(cmap, range, center):
    """
    Converts center to a ratio between 0 and 1 of the
    range given and calls cmap_center_adjust(). returns
    a new adjusted colormap accordingly

    :param cmap: colormap instance
    :param range: Tuple of (min, max)
    :param center: New cmap center
    """
    if not ((range[0] < center) and (center < range[1])):
        return cmap
    return cmap_center_adjust(
        cmap,
        abs(center - range[0]) / abs(range[1] - range[0]))


if __name__ == '__main__':
    def func3(x, y):
        return (1 - x / 2 + x ** 5 + y ** 3) * numpy.exp(-x ** 2 - y ** 2)
    x = numpy.linspace(-3.0, 3.0, 60)
    y = numpy.linspace(-3.0, 3.0, 60)
    X, Y = numpy.meshgrid(x, y)
    Z = func3(X, Y)
    extent = [x[0], x[-1], y[0], y[-1]]

    plotkwargs = {
        'extent': extent,
        'origin': 'lower',
        'interpolation': 'nearest',
        'aspect': 'auto'}

    fig = pyplot.figure(figsize=(8, 3))
    fig.subplots_adjust(left=.05, bottom=.11, right=.94, top=.83, wspace=.35)
    ax = [fig.add_subplot(1, 3, i) for i in range(1, 4, 1)]

    cmap = cm.seismic

    plt = ax[0].imshow(Z, cmap=cmap, **plotkwargs)
    cb = ax[0].figure.colorbar(plt, ax=ax[0])
    ax[0].set_title('cmap: seismic')

    plt = ax[1].imshow(Z, cmap=cmap_center_adjust(cmap, 0.75), **plotkwargs)
    cb = ax[1].figure.colorbar(plt, ax=ax[1])
    ax[1].set_title('center raised by 25%')

    plt = ax[2].imshow(
        Z,
        cmap=cmap_center_point_adjust(
            cmap, [numpy.min(Z), numpy.max(Z)], 0),
        **plotkwargs)
    cb = ax[2].figure.colorbar(plt, ax=ax[2])
    ax[2].set_title('center set to 0')

    pyplot.show()
