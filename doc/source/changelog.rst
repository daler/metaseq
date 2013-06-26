Changelog
=========
Changes in 0.5dev
-----------------
* rebinning uses interpolation, rather than trying to keep everything as
  integers, resulting in far fewer artifacts when doing genome-wide averages

* better bigwig support

* `results_table` object now wraps `pandas.DataFrame`, so is much more flexible

* lots of bugfixes, lots and lots of API improvements
