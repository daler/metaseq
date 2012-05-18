API docs
========
.. contents::


.. currentmodule:: metaseq.genomic_signal
:mod:`genomic_signal`
---------------------
.. automodule:: metaseq.genomic_signal

----

.. autosummary::
   :nosignatures:

   metaseq.genomic_signal.genomic_signal
   metaseq.genomic_signal.supported_formats
   metaseq.genomic_signal.BigWigSignal
   metaseq.genomic_signal.BamSignal
   metaseq.genomic_signal.BigBedSignal
   metaseq.genomic_signal.BedSignal


----

.. autofunction:: metaseq.genomic_signal.genomic_signal

----

.. autofunction:: metaseq.genomic_signal.supported_formats

----

.. autoclass:: BigWigSignal
    :show-inheritance:
    :inherited-members:
    :members:

----

.. autoclass:: BamSignal
    :show-inheritance:
    :inherited-members:
    :members:

----

.. autoclass:: BigBedSignal
    :show-inheritance:
    :inherited-members:
    :members:

----

.. autoclass:: BedSignal
    :show-inheritance:
    :inherited-members:
    :members:


:mod:`integration`
------------------
.. currentmodule:: metaseq.integration

Module that ties together the various parts of :mod:`metaseq`.

----

.. autosummary::
    :nosignatures:

    metaseq.integration.chipseq.Chipseq

:mod:`integration.chipseq`
~~~~~~~~~~~~~~~~~~~~~~~~~~
.. automodule:: metaseq.integration.chipseq

----

.. autoclass:: metaseq.integration.chipseq.Chipseq




:mod:`metaseq.minibrowser`
--------------------------
.. currentmodule:: metaseq.minibrowser


.. autosummary::
    :nosignatures:

    metaseq.minibrowser.BaseMiniBrowser
    metaseq.minibrowser.SignalMiniBrowser
    metaseq.minibrowser.GeneModelMiniBrowser

.. automodule:: metaseq.minibrowser

----

.. autoclass:: metaseq.minibrowser.BaseMiniBrowser
    :show-inheritance:
    :members:

----

.. autoclass:: metaseq.minibrowser.SignalMiniBrowser
    :show-inheritance:
    :members:

----

.. autoclass:: metaseq.minibrowser.GeneModelMiniBrowser
    :show-inheritance:
    :members:


:mod:`metaseq.filetype_adapters`
--------------------------------
.. automodule:: metaseq.filetype_adapters
.. autosummary::
    :nosignatures:

    metaseq.filetype_adapters.BaseAdapter
    metaseq.filetype_adapters.BamAdapter
    metaseq.filetype_adapters.BedAdapter
    metaseq.filetype_adapters.BigBedAdapter
