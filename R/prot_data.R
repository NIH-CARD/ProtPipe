#' ProtData Class
#'
#' An S4 class that holds proteomics data and provides methods for processing.
#'
#' @export
setClass("ProtData",
         slots = list(
           data = "data.frame",      # The main proteomics data (proteins are rows, samples are columns)
           condition = "data.frame",  # Conditions of the samples. Rownames must match colnames of data
           prot_meta = "data.frame",  # information about the rows (genes, organism, etc...)
           method = "character"       # The method used (e.g., "MS", "Somascan", etc.)
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
create_protdata <- function(data, condition = NULL, prot_meta = NULL, method = "Unknown") {

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

  # Create a new ProtData object
  new("ProtData",
      data = data,
      condition = condition,
      prot_meta = prot_meta,
      method = method)
}
