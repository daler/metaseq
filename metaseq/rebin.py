import numpy as np
def rebin(x, y, nbin):
    xi = np.linspace(x.min(), x.max(), nbin)
    return xi, np.interp(xi, x, y)
