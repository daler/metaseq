Changelog
=========
Changes in v0.5.5.1
-------------------
* change DifferentialExpressionResults.scatter default plotting style to have
  linewidths > 0.  Previously, using linewidths=0 cause plots to show up blank
  in IPython Notebook on Mac

Changes in v0.5.5
-----------------
* add `metaseq.tableprinter` module

Changes in v0.5.4
-----------------
* `metaseq.results_table.DESeq2Results` class
* handle cases where callbacks could be called multiple times
* support marginal histograms of length-1 arrays
* `metaseq.persistence` functions use absolute path for features to make
  symlinking more robust
* differential expression classes cleaned up and refactored
* new kwargs for `metaseq.plotutils.imshow`: `figure_kwargs`, `subset_by`,
  `subset_order`
* new `metaseq.plotutils.add_labels_to_subsets` helper function


Changes in v0.5
---------------
v0.5 marks a major improvement and maturation of `metaseq` but at the expense
of backward-compatibility.  The changes are too numerous to list, but include
massive plotting improvements, speedups, more flexible API, and much more.
