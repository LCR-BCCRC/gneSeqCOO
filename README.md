# gneSeqCOO

This package contains code for applying a Cell of Origin classifier to RNASeq data. The process consists of three steps: 
	1) normalizing
	2) calculating the LPS score
	3) splitting samples into GCB/ABC based on fixed cut-points.

# Install

```
remotes::install_github("Genentech/gneSeqCOO")
```

# Basic usage

The input for the gneSeqCOO method is a raw count matrix, formatted as a DESeqDataSet. For a standard count matrix, this can be generated via the code below:

```{r,eval=FALSE}
## counts is a matrix, with rows of features and columns of samples.
cds = DESeqDataSetFromMatrix(counts,
  colData=data.frame(ID=colnames(counts)),
  rowData=data.frame(ID=rownames(counts)),
  design=~1)
```

Note that in order for the classifier to work, it is necessary for the row names of our CDS to be in the correct format. The classifier currently accepts Refseq “Gene IDs” (formatted as “GeneID:9294”) or Ensembl IDs (ENSG....). All genes used in the model must be present in the dataset. For the standard model, there are 21 genes, with the following IDs:

```{r,echo=FALSE,message=FALSE}
library(gneSeqCOO)
themodel=coo_GOYA_21gene_igis3$coefficients
names(themodel)[names(themodel)!="(Intercept)"]
```

The COO classifier can be applied to the cds object using the single command, coo_rnaseq(). An example is below:

```{r,eval=FALSE}

pred = coo_rnaseq(cds)

```

The returned object consists of three columns:
  - Sample IDs, drawn from the column names of the cds object.
  - LPS - the Linear Predictor Score, used to split samples into GCB and ABC, and
  - COO - the Cell of Origin classification; either GCB, ABC, or Unclassified


For more details, please refer to the docmumentation in the package.
