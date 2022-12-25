# Summary

We used a DNA micro-array chip covering 369 resistance types to investigate the relation of antibiotic resistance gene diversity with humans' age.

**Meta-genomic DNA from fecal samples of 123 healthy volunteers of four different age groups:**

-   pre-school Children (CH)

-   School Children (SC)

-   High School Students (HSS)

-   Adults (AD)

**The results showed that**

-   80 different gene types were recovered from the 123 individuals gut microbiota, among which 25 were present in CH, 37 in SC, 58 in HSS and 72 in AD.

-   Further analysis indicated that antibiotic resistance genes in groups of CH, SC and AD can be independently clustered, and those ones in group HSS are more divergent.

------------------------------------------------------------------------

## Code clarification:

-   **Normalization**
    -   Normalization of micro-array data typically consists of a series of data transformations intended to aid in the comparison of gene expression data gathered across a series of hybridization.

    -   These may include background correction to remove geographical biases in fluorescent intensity

        -   **applying intensity thresholds or flooring** to remove poorly detected probes and to increase signal-noise sensitivity

        -   **log transformation** to normalize the distribution of probes across the intensity range of the experiment, and the scaling of data such that information extracted from one array slide is equivalent to that extracted from another in the series.

    -   The process of threshold, scaling and log transforming data will reduce variance between samples.
-   **Box-and-whisker plot**
    -   **par() function =\>** Combining multiple graphs into a single image for our viewing convenience.

    -   **boxwex =\>** a scale factor to be applied to all boxes. When there are only a few groups, the appearance of the plot can be improved by making the boxes narrower.

    -   **las =\>** Orient the axis labels for bplot groups. Default is to put them horizontal (las=1) if number of groups is less than 7 otherwise make them perpendicular (las=2) to keep the labels from running into each other. See help(par) for more about this option.
-   **Expression value distribution**
-   **Variance of the mean**
    -   Generates a subsample variance of the mean versus subsample index plot.

    -   The subsample variance of the mean is the subsample variance divided by the subsample size.

    -   The variance of the mean plot is used to answer the question: "Does the subsample variation of the mean change over different subsamples?"

    -   **plotSA**

        -   Plot residual standard deviation versus average log expression for a fitted microarray linear model.

        -   calculates how much the data points spread around the regression line.
