################################################################################
# Updated version of the code for the analysis in:
#
#   "Multivariate meta-analysis for non-linear and other 
#     multi-parameter associations"
#   Gasparrini, Armstrong and Kenward
#   Statistics in Medicine 2012
#   http://www.ag-myresearch.com/statmed2012.html
#
# Update: 11 April 2016
# For any problem with this code, please contact antonio.gasparrini@lshtm.ac.uk
# Please refer to the original code for any copyright issue
#
#  See www.ag-myresearch.com for future updates
################################################################################

The original example included in the article was based on data from the National Mortality, Morbidity, and Air Pollution Study (NMMAPS), which at the time of the publication was available through the R package NMMAPSlite.

Unfortunately, the data are not available any more and the package NMMAPSlite has been archived. This means that the analysis of the paper is not replicable.

In order to provide a working example, the code has been replaced with a similar analysis on a publicly available dataset, with series of mortality for 10 regions in England and Wales for the period 1993-2006. The data are stored in the file regEngWales.csv. The same data are also used in another article available at http://www.ag-myresearch.com/bmcmrm2013.html, where additional information are also available.

A bug was fixed on 22 February 2016. Thanks to  Felipe J Colon-Gonzalez.
