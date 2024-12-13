#' ProtData Class
#'
#' An S4 class that holds proteomics data and provides methods for processing.
#'
#' @export
setClass("ProtData",
         slots = list(
           data = "data.table",      # The main proteomics data (proteins are rows, samples are columns)
           data.long = "data.table",
           condition = "data.table",  # Conditions of the samples. Rownames must match colnames of data
           prot_meta = "data.table",  # information about the rows (genes, organism, etc...)
           method = "character"       # The method used (e.g., "MS", "Somascan", etc.)
         ),
         prototype = list(
           data = data.table::data.table(),        # Default is an empty data frame
           data.long = data.table::data.table(),   # Default is an empty data frame
           condition = data.table::data.table(),   # Default is an empty data frame
           prot_meta = data.table::data.table(),   # Default is an empty data frame
           method = "Unknown"          # Default method is set to "Unknown"
         )
)



# Constructor for ProtData class
#' Create a ProtData Object
#'
#' This function creates an instance of the ProtData class.
#'
#' @param data A data frame containing proteomics data (proteins are rows, samples are columns).
#' @param condition A data frame containing conditions of the samples. Rownames should match colnames of data. Optional.
#' @param prot_meta A data frame containing metadata for the proteins (rows). Optional.
#' @param method A character string describing the method used for generating the data. Optional.
#'
#' @return An instance of the ProtData class.
#' @export
create_protdata <- function(data, condition = data.table::data.table(), prot_meta = data.table::data.table(), method = "Unknown") {

  # Check that data is a data frame
  if (!is.data.frame(data)) {
    stop("The 'data' argument must be a data frame.")
  }

  # Ensure that condition, if provided, has rownames matching the colnames of data
  if (!is.null(condition)) {
    if (!all(rownames(condition) %in% colnames(data))) {
      stop("Rownames of 'condition' must match the colnames of 'data'.")
    }
  }
  #standardize the column names and order
  data <- standardize_format(data)
  data.table::setnames(dat, trim_colnames(dat))
  col_order=c(colnames(dat)[1:2],sort(colnames(dat)[3:ncol(dat)]))
  data.table::setcolorder(dat,col_order)
  data <- data %>%
    # Calculate missing values
    dplyr::mutate(missing_value = rowSums(is.na(dplyr::select(., -contains("Protein_Group|Genes")))))%>%
    # Calculate median values
    dplyr::mutate(median = matrixStats::rowMedians(as.matrix(dplyr::select(., -contains("Protein_Group|Genes|missing_value")) %>%
                                           dplyr::select_if(is.numeric)), na.rm = TRUE)) %>%
    #Identify unique Protein_Group with the least missing values and highest median intensity
    dplyr::group_by(Protein_Group) %>%
    dplyr::filter(missing_value == min(missing_value)) %>%
    dplyr::slice(which.max(median)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-c(missing_value, median))

  dat.long <- melt_intensity_table(data)
  dat.long <- dat.long[rowSums(is.na(dat.long)) < ncol(dat.long),]
  #dat.long <- dat.long[! is.na(Intensity)][Intensity != 0]
  #filter intensity
  #dat.long <- dat.long[Intensity > opt$minintensity]

  # Create a new ProtData object
  new("ProtData",
      data = data.table::setDT(data),
      data.long = dat.long,
      condition = condition,
      prot_meta = prot_meta,
      method = method)
}


####### GETTERS and SETTERS ###############################################################

# Define the generic for 'getData'
setGeneric("getData", function(object) standardGeneric("getData"))

# Define the setter for 'getData'
setGeneric("setData", function(object, value) standardGeneric("setData"))

# Now define the method for 'getData' for the ProtData class
setMethod("getData", "ProtData", function(object) {
  return(object@data)
})

# Define the setter method for 'setData' for ProtData
setMethod("setData", "ProtData", function(object, value) {
  if (!inherits(value, "data.table")) {
    stop("'data' must be a data.table.")
  }
  object@data <- value
  return(object)
})

# Define the generic for 'data.long'
setGeneric("getDataLong", function(object) standardGeneric("getDataLong"))
setGeneric("setDataLong", function(object, value) standardGeneric("setDataLong"))

# Define the methods for 'data.long' for ProtData
setMethod("getDataLong", "ProtData", function(object) {
  return(object@data.long)
})

setMethod("setDataLong", "ProtData", function(object, value) {
  if (!inherits(value, "data.table")) {
    stop("'data.long' must be a data.table.")
  }
  object@data.long <- value
  return(object)
})

# Similarly, define getter and setter for 'condition'
setGeneric("getCondition", function(object) standardGeneric("getCondition"))
setGeneric("setCondition", function(object, value) standardGeneric("setCondition"))

setMethod("getCondition", "ProtData", function(object) {
  return(object@condition)
})

setMethod("setCondition", "ProtData", function(object, value) {
  if (!inherits(value, "data.table")) {
    stop("'condition' must be a data.table.")
  }
  object@condition <- value
  return(object)
})

# Similarly, define getter and setter for 'prot_meta'
setGeneric("getProtMeta", function(object) standardGeneric("getProtMeta"))
setGeneric("setProtMeta", function(object, value) standardGeneric("setProtMeta"))

setMethod("getProtMeta", "ProtData", function(object) {
  return(object@prot_meta)
})

setMethod("setProtMeta", "ProtData", function(object, value) {
  if (!inherits(value, "data.table")) {
    stop("'prot_meta' must be a data.table.")
  }
  object@prot_meta <- value
  return(object)
})

# Similarly, define getter and setter for 'method'
setGeneric("getProtMethod", function(object) standardGeneric("getProtMethod"))
setGeneric("setProtMethod", function(object, value) standardGeneric("setProtMethod"))

#WHY DOESNT THIS WORK????
setMethod("getProtMethod", "ProtData", function(object) {
  return(object@method)
})

setMethod("setProtMethod",
          signature = c(object = "ProtData", value = "character"),
          function(object, value) {
            if (length(value) != 1) {
              stop("'method' must be a single character string.")
            }
            object@method <- value
            return(object)
          })







######## HELPER METHODS ####################################################################

#total proteomics############
standardize_format <- function(DT.original) {
  # Accepts an input protein group intensity data.table, whether spectronaut or DIA-NN format,
  # and restructures into one consistent style for downstream processing
  DT <- DT.original
  if("Protein.Ids" %in% colnames(DT)) {
    print("DIAnn input")
    DT[, 'Protein.Ids' := NULL]
    DT[, 'Protein.Names' := NULL]
    DT[, 'First.Protein.Description' := NULL]
    data.table::setnames(DT, 'Protein.Group', 'Protein_Group')
  }
  else if('EG.PrecursorId' %in% colnames(DT)) {
    print("Spectronaut input")
    data.table::setnames(DT, 'EG.PrecursorId', 'Peptide_Sequence')
    data.table::setnames(DT, 'PG.Genes', 'Genes')
    DT=as.data.frame(DT)
    # Use only Protein_Group and Genes
    col_dplyr::select=c('Peptide_Sequence','Genes',grep('raw',colnames(DT),value = T))
    DT=DT[, col_dplyr::select]
    #as number
    DT[,grep('raw',colnames(DT))]=as.data.frame(apply(DT[,grep('raw',colnames(DT))],2,as.numeric))
    DT=data.table(DT)
  }
  else if('PG.ProteinGroups' %in% colnames(DT)) {
    print("Spectronaut input")
    data.table::setnames(DT, 'PG.ProteinGroups', 'Protein_Group')
    data.table::setnames(DT, 'PG.Genes', 'Genes')
    DT=as.data.frame(DT)
    # Use only Protein_Group and Genes
    col_dplyr::select=c('Protein_Group','Genes',grep('PG.Quantity',colnames(DT),value = T))
    DT=DT[, col_dplyr::select]
    #as number
    DT[,grep('PG.Quantity',colnames(DT))]=as.data.frame(apply(DT[,grep('PG.Quantity',colnames(DT))],2,as.numeric))
    DT=data.table(DT)
  }
  else if('Peptide Sequence' %in% colnames(DT)) {
    print("FragPipe input")
    data.table::setnames(DT, 'Gene', 'Genes')
    colnames(DT)=gsub("\\s", "_",colnames(DT))
    # Use only Protein_Group and Genes
    DT=as.data.frame(DT)
    col_dplyr::select=c('Peptide_Sequence','Genes',grep('[0-9]_Intensity',colnames(DT),value = T))
    DT=DT[, col_dplyr::select]
    DT=data.table(DT)
  }

  # Remove leading directories for sample names
  # e.g. /path/to/sample1.mzML -> sample1.mzML
  data.table::setnames(DT, basename(colnames(DT)))

  # Remove trailing file extensions
  extensions <- '.mzML$|.mzml$|.RAW$|.raw$|.dia$|.DIA$|_Intensity'
  extension_samplenames <-  colnames(DT)[data.table::`%like%`(colnames(DT), extensions)]
  trimmed_samplenames <- gsub(extensions, '', extension_samplenames)
  data.table::setnames(DT, extension_samplenames, trimmed_samplenames)
  return(DT[])
}

trim_colnames <- function(DT) {
  colnames_out <- gsub(pattern="\\[.*\\] ", replacement='', x=colnames(DT))   # trim leading [N]
  colnames_out <- gsub(pattern="\\..*\\.PG\\.Quantity|\\.PG\\.Quantity|\\..*Quantity.*", replacement='', x=colnames_out)   # remove suffix
  return(colnames_out)
}

#create long data table
melt_intensity_table <- function(DT) {
  # Converts intensity data.table to long format
  # info_cols <- c('Protein_Group', 'Genes', 'First_Protein_Description')
  DT.long <- reshape2::melt(DT,
                  measure.vars=names(DT)[sapply(DT, function(x) all(is.numeric(x)))],
                  variable.name='Sample',
                  value.name='Intensity')
  DT.long=data.table::data.table(DT.long)
  DT.long=DT.long[DT.long$Intensity>0,]
  return(DT.long)
}
