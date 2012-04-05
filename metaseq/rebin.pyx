import numpy as np
cimport numpy as np
cimport cython
DTYPE_INT = np.int
DTYPE_FLOAT = np.float
ctypedef np.int_t DTYPE_INT_t
ctypedef np.float_t DTYPE_FLOAT_t

"""
Optimization notes:

    This version is 6-10x faster than the naive Python version.

    Counterintuitively, np.sum was slow, so the sums were re-implemented using
    for-loops iterating over indexes.  Since that's rewritten to C loops, it's
    faster.

    y = np.array([0,  0,  1,  2,  2,  1,  1,  0,  0,  0])

rebin.rebin
-----------

    timeit rebin.rebin(y, 10)
    100000 loops, best of 3: 2.96 us per loop

    timeit rebin.rebin(y, 20)
    100000 loops, best of 3: 8.1 us per loop

    timeit rebin.rebin(y, 5)
    100000 loops, best of 3: 3 us per loop


naive_rebin
-----------

    timeit naive_rebin(y, 5)
    10000 loops, best of 3: 30.4 us per loop

    timeit naive_rebin(y, 10)
    10000 loops, best of 3: 51.3 us per loop

    timeit f(y, 20)
    10000 loops, best of 3: 23.5 us per loop

signal.resample
---------------

    timeit signal.resample(y, 5)
    10000 loops, best of 3: 55.4 us per loop

    timeit signal.resample(y, 10)
    10000 loops, best of 3: 56.3 us per loop

    timeit signal.resample(y, 20)
    10000 loops, best of 3: 56.8 us per loop


use a larger array
-----------------

    timeit f(arr, 100)
    1000 loops, best of 3: 506 us per loop

    timeit rebin.rebin(arr, 100)
    10000 loops, best of 3: 33.4 us per loop

    timeit signal.resample(arr, 100)
    1000 loops, best of 3: 251 us per loop


"""
def float_rebin(np.ndarray[DTYPE_FLOAT_t, ndim=1] y, int bins):
    """
    Given a numpy array of integers, `y`, re-distribute the items into `bins`
    equally-spaced bins.
    """
    cdef np.ndarray[DTYPE_FLOAT_t, ndim=1] yi = np.zeros( (bins,), dtype=DTYPE_FLOAT)
    cdef int i, j, i1, i2, padding, c
    cdef DTYPE_INT_t s
    cdef int len_y = y.shape[0]

    # exit quickly if there's nothing to do
    if len_y == bins:
        return y

    # this determines how many indices in `y` should be summed into a single
    # bin
    cdef int spacing  = int(np.round(len_y / float(bins)))

    if spacing > 0:
        for i in xrange(bins):
            i1 = i * spacing
            i2 = i1 + spacing
            for k in xrange(i1, i2):
                # make sure we don't overshoot y
                if k < len_y:
                    yi[i] += y[k]

        # If anything in `y` hasn't been added yet (from rounding in the
        # `spacing` calculation) then add the remainder to the last bin
        if i2 < len_y:
            for k in xrange(i2, len_y - 1):
                yi[bins - 1] += y[k]

    # If spacing == 0, then we need to expand `y` rather than collapse it into
    # fewer bins.  A choice to make here:  
    # 1) leaving the expanded bins == 0 means the signal will not be continuous,
    #    but it guarantees sum(yi) == sum(y)
    # 2) filling in the intervening bins means the signal will look smooth, but
    #    sum(yi) > sum(yi).
    else:
        # This is how many indicies to put after each actual value in `y`.
        # Another way of looking at it is that it's the number of zeros to
        # append to `yi` after appending each value in `y`
        padding = int(np.floor(1. / (float(len_y) / float(bins))))

        c = 0
        for i in xrange(len_y):
            yi[c] = y[i]
            c += 1
            for j in xrange(padding - 1):
                if c < bins:
                    # Uncomment to use option 2 above
                    #yi[c] = y[i]
                    c += 1

    return yi

def rebin(np.ndarray[DTYPE_INT_t, ndim=1] y, int bins):
    """
    Given a numpy array of integers, `y`, re-distribute the items into `bins`
    equally-spaced bins.
    """
    cdef np.ndarray[DTYPE_INT_t, ndim=1] yi = np.zeros( (bins,), dtype=DTYPE_INT)
    cdef int i, j, i1, i2, padding, c
    cdef DTYPE_INT_t s
    cdef int len_y = y.shape[0]

    # exit quickly if there's nothing to do
    if len_y == bins:
        return y

    # this determines how many indices in `y` should be summed into a single
    # bin
    cdef int spacing  = int(np.round(len_y / float(bins)))

    if spacing > 0:
        for i in xrange(bins):
            i1 = i * spacing
            i2 = i1 + spacing
            for k in xrange(i1, i2):
                # make sure we don't overshoot y
                if k < len_y:
                    yi[i] += y[k]

        # If anything in `y` hasn't been added yet (from rounding in the
        # `spacing` calculation) then add the remainder to the last bin
        if i2 < len_y:
            for k in xrange(i2, len_y - 1):
                yi[bins - 1] += y[k]

    # If spacing == 0, then we need to expand `y` rather than collapse it into
    # fewer bins.  A choice to make here:  
    # 1) leaving the expanded bins == 0 means the signal will not be continuous,
    #    but it guarantees sum(yi) == sum(y)
    # 2) filling in the intervening bins means the signal will look smooth, but
    #    sum(yi) > sum(yi).
    else:
        # This is how many indicies to put after each actual value in `y`.
        # Another way of looking at it is that it's the number of zeros to
        # append to `yi` after appending each value in `y`
        padding = int(np.floor(1. / (float(len_y) / float(bins))))

        c = 0
        for i in xrange(len_y):
            yi[c] = y[i]
            c += 1
            for j in xrange(padding - 1):
                if c < bins:
                    # Uncomment to use option 2 above
                    #yi[c] = y[i]
                    c += 1

    return yi
