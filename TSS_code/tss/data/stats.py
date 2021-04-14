from scipy.interpolate import splrep, sproot, splev
import numpy as np


def fwhm(x, y, k=3):
    """
    Determine full-width-half-maximum of a peaked set of points,
    x and y.

    Assumes that there is only one peak present in the datasset.  The function
    uses a spline interpolation of order k.

    Reference: https://stackoverflow.com/questions/10582795/finding
    -the-full-width-half-maximum-of-a-peak/16489955 User jdg
    """

    class MultiplePeaks(Exception):
        pass

    class NoPeaksFound(Exception):
        pass

    half_max = np.max(y)/2.0

    print(half_max)
    s = splrep(x, y - half_max, k=k)
    roots = sproot(s)
    print(roots)
    if len(roots) > 2:
        raise MultiplePeaks("The dataset appears to have multiple peaks, and "
                "thus the FWHM can't be determined.")
    elif len(roots) < 2:
        raise NoPeaksFound("No proper peaks were found in the data set; likely "
                "the dataset is flat (e.g. all zeros).")
    else:
        return abs(roots[1] - roots[0])

