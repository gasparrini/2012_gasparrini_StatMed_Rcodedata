### Multivariate meta-analysis for non-linear and other multi-parameter associations

------------------------------------------------------------------------

The development and application of multivariate meta-analysis for pooling estimates of non-linear associations from multiple studies. The code illustrates an example in the specific setting of time series analysis of temperature-health relationships, but the methodology is generally applicable in a broader context. The example illustrates the analysis in the article:

Gasparrini A, Armstrong B, Kenward MG. Multivariate meta-analysis for non-linear and other multi-parameter associations. *Statistics in Medicine*. 2012;**31**(29):3821-39. DOI: 10.1002/sim.5471. PMID: 22807043. [[freely available here](http://www.ag-myresearch.com/2012_gasparrini_statmed.html)]

The original example included in the article was based on data from the National Mortality, Morbidity, and Air Pollution Study (NMMAPS), which at the time of the publication was available through the R package NMMAPSlite. Unfortunately, the data are not available any more and the package NMMAPSlite has been archived. In order to provide a working example, the analysis has been replaced with a similar analysis on a publicly available dataset.

The code primarily uses functions in the [R package mvmeta](https://github.com/gasparrini/mvmeta), but functions in the [R package dlnm](https://github.com/gasparrini/dlnm) are applied as well.

------------------------------------------------------------------------

The material:

-   *regEngWales.csv* stores the dataset used to perform an illustrative example, including time series data for 10 regions in England and Wales for 1993-2006
-   the numbered files from *01.mkdata.R* to *06.addresults.R*, to be run sequentially, reproduce the example (different from that in the article)

Download as a ZIP file using the green button *Clone or download* above
