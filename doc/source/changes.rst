Changelog
=========
Changes in v0.5.6
-----------------
(versions 0.5.5.2-4 were minor versions involving different testing frameworks)

* fix typo in plotutils.clustered_sortind
* support for versions of pysam expecting string (rather than unicode) for
  chromosome
* support for cells in 2x2 tables with 0 items
* lots of improvements to testing framework -- now using docker on travis
* added fisher as an explicit requirement
* handle bx-python 0.7.3 which does not work correctly on bigbed files.
  A NotImplementedError is raised in this case. The workaround is to convert
  bigBed to BED (`bigBedToBed` program from UCSC; `conda install -c bioconda
  ucsc-bigbedtobed` to install) and then either use a tabix-indexed BED file or
  convert BED to BAM (`bedtools bedtobam` program).


Changes in v0.5.5.1
-------------------
* change DifferentialExpressionResults.scatter default plotting style to have
  linewidths > 0.  Previously, using linewidths=0 caused plots to show up blank
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
