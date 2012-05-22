API docs
========
.. contents::


.. currentmodule:: metaseq.genomic_signal

:mod:`genomic_signal`
---------------------
.. automodule:: metaseq.genomic_signal

----

.. rubric:: Functions:

.. autosummary::
   :nosignatures:
   :toctree: autodocs

   metaseq.genomic_signal.genomic_signal
   metaseq.genomic_signal.supported_formats

.. rubric:: Classes

.. autosummary::
   :nosignatures:
   :toctree: autodocs
   :template: auto_template.rst

   metaseq.genomic_signal.BigWigSignal
   metaseq.genomic_signal.BamSignal
   metaseq.genomic_signal.BigBedSignal
   metaseq.genomic_signal.BedSignal


----

:mod:`metaseq.integration`
--------------------------

.. automodule:: metaseq.integration

----

.. rubric:: Classes

.. autosummary::
    :nosignatures:
    :toctree: autodocs
    :template: auto_template.rst

    metaseq.integration.chipseq.Chipseq

.. rubric:: Functions

.. autosummary::
    :nosignatures:
    :toctree: autodocs

    metaseq.integration.signal_comparison.compare

:mod:`metaseq.colormap_adjust`
------------------------------
.. automodule:: metaseq.colormap_adjust

.. autosummary::
    :nosignatures:
    :toctree: autodocs

    metaseq.colormap_adjust.color_test
    metaseq.colormap_adjust.smart_colormap
    metaseq.colormap_adjust.cmap_discretize
    metaseq.colormap_adjust.cmap_powerlaw_adjust
    metaseq.colormap_adjust.cmap_center_adjust
    metaseq.colormap_adjust.cmap_center_point_adjust

----

:mod:`metaseq.minibrowser`
--------------------------
.. currentmodule:: metaseq.minibrowser

.. automodule:: metaseq.minibrowser

.. autosummary::
    :nosignatures:
    :toctree: autodocs
    :template: auto_template.rst

    metaseq.minibrowser.BaseMiniBrowser
    metaseq.minibrowser.SignalMiniBrowser
    metaseq.minibrowser.GeneModelMiniBrowser





:mod:`metaseq.filetype_adapters`
--------------------------------
.. automodule:: metaseq.filetype_adapters
.. autosummary::
    :nosignatures:
    :toctree: autodocs
    :template: auto_template.rst

    metaseq.filetype_adapters.BaseAdapter
    metaseq.filetype_adapters.BamAdapter
    metaseq.filetype_adapters.BedAdapter
    metaseq.filetype_adapters.BigBedAdapter
