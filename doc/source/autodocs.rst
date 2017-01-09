API docs
========
.. contents::



:mod:`genomic_signal`
---------------------
.. automodule:: metaseq._genomic_signal

----

.. rubric:: Functions:

.. autosummary::
   :nosignatures:
   :toctree: autodocs

   metaseq._genomic_signal.genomic_signal
   metaseq._genomic_signal.supported_formats

.. rubric:: Classes

.. autosummary::
   :nosignatures:
   :toctree: autodocs
   :template: auto_template.rst

   metaseq._genomic_signal.BaseSignal
   metaseq._genomic_signal.IntervalSignal
   metaseq._genomic_signal.BigWigSignal
   metaseq._genomic_signal.BamSignal
   metaseq._genomic_signal.BigBedSignal
   metaseq._genomic_signal.BedSignal


----

:mod:`metaseq.results_table`
----------------------------
.. automodule:: metaseq.results_table

.. rubric:: Classes


.. autosummary::
   :nosignatures:
   :toctree: autodocs
   :template: auto_template.rst

   metaseq.results_table.ResultsTable
   metaseq.results_table.DESeqResults
   metaseq.results_table.DESeq2Results
   metaseq.results_table.EdgeRResults
   metaseq.results_table.LazyDict

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

----

:mod:`metaseq.plotutils`
------------------------

.. automodule:: metaseq.plotutils

.. rubric:: Functions

.. autosummary::
    :nosignatures:
    :toctree: autodocs

    metaseq.plotutils.imshow
    metaseq.plotutils.add_labels_to_subsets
    metaseq.plotutils.clustered_sortind
    metaseq.plotutils.calculate_limits
    metaseq.plotutils.ci_plot
    metaseq.plotutils.ci
    metaseq.plotutils.tip_zscores
    metaseq.plotutils.tip_fdr
    metaseq.plotutils.nice_log
    metaseq.plotutils.prepare_logged
    metaseq.plotutils.matrix_and_line_shell
    metaseq.plotutils.input_ip_plots

.. rubric:: Classes

.. autosummary::
    :nosignatures:
    :toctree: autodocs
    :template: auto_template.rst

    metaseq.plotutils.MarginalHistScatter

----

:mod:`metaseq.colormap_adjust`
------------------------------
.. automodule:: metaseq.colormap_adjust

.. rubric:: Functions

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

.. automodule:: metaseq.minibrowser

.. rubric:: Classes

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

.. rubric:: Classes

.. autosummary::
    :nosignatures:
    :toctree: autodocs
    :template: auto_template.rst

    metaseq.filetype_adapters.BaseAdapter
    metaseq.filetype_adapters.BamAdapter
    metaseq.filetype_adapters.BedAdapter
    metaseq.filetype_adapters.BigBedAdapter
    metaseq.filetype_adapters.BigWigAdapter

:mod:`metaseq.array_helpers`
----------------------------
.. automodule:: metaseq.array_helpers

.. rubric:: Functions

.. autosummary::
    :nosignatures:
    :toctree: autodocs

    _local_count
    _local_coverage

