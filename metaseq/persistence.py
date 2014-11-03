import os
import pybedtools
import numpy as np

"""
Tools for working with data across sessions.
"""


def load_features_and_arrays(prefix, mmap_mode='r'):
    """
    Returns the features and NumPy arrays that were saved with
    save_features_and_arrays.

    Parameters
    ----------

    prefix : str
        Path to where data are saved

    mmap_mode : {None, 'r+', 'r', 'w+', 'c'}
        Mode in which to memory-map the file.  See np.load for details.
    """
    features = pybedtools.BedTool(prefix + '.features')
    arrays = np.load(prefix + '.npz', mmap_mode=mmap_mode)
    return features, arrays


def save_features_and_arrays(features, arrays, prefix, compressed=False,
                             link_features=False, overwrite=False):
    """
    Saves NumPy arrays of processed data, along with the features that
    correspond to each row, to files for later use.

    Two files will be saved, both starting with `prefix`:

        prefix.features : a file of features.  If GFF features were provided,
        this will be in GFF format, if BED features were provided it will be in
        BED format, and so on.

        prefix.npz : A NumPy .npz file.

    Parameters
    ----------
    arrays : dict of NumPy arrays
        Rows in each array should correspond to `features`.  This dictionary is
        passed to np.savez

    features : iterable of Feature-like objects
        This is usually the same features that were used to create the array in
        the first place.

    link_features : bool
        If True, then assume that `features` is either a pybedtools.BedTool
        pointing to a file, or a filename.  In this case, instead of making
        a copy, a symlink will be created to the original features.  This helps
        save disk space.

    prefix : str
        Path to where data will be saved.

    compressed : bool
        If True, saves arrays using np.savez_compressed rather than np.savez.
        This will save disk space, but will be slower when accessing the data
        later.
    """

    if link_features:
        if isinstance(features, pybedtools.BedTool):
            assert isinstance(features.fn, basestring)
            features_filename = features.fn
        else:
            assert isinstance(features, basestring)
            features_filename = features

        if overwrite:
            force_flag = '-f'
        else:
            force_flag = ''

        cmds = [
            'ln', '-s', force_flag, os.path.abspath(features_filename), prefix + '.features']
        os.system(' '.join(cmds))
    else:
        pybedtools.BedTool(features).saveas(prefix + '.features')

    if compressed:
        np.savez_compressed(
            prefix,
            **arrays)
    else:
        np.savez(prefix, **arrays)
