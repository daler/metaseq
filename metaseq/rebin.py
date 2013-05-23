import numpy as np
def rebin(y, bins):
    len_y = y.shape[0]
    x = np.arange(len_y)
    xi = np.linspace(0, len_y, bins)
    return np.interp(xi, x, y)
