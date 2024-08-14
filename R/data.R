#' Primary COO classifier
#'
#' An 'lm' object containing a 21 gene COO classifier, developed from the GOYA 
#' dataset. Two versions are available, for data processed with either
#' IGIS3.0 or IGIS4.0
#'
#' @name coo_GOYA_21gene
#' @format An object of class 'lm'. Coefficients can be viewed by accessing
#'    the '$coefficients' slot.
#' @source Refer to the COO_Classifier.html file in inst/reports for details
#' on the development of this model.
"coo_GOYA_21gene_igis3"


#' @rdname coo_GOYA_21gene
"coo_GOYA_21gene_igis4"