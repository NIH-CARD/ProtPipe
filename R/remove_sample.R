#' Title
#'
#' @param object
#' @param samples
#'
#' @return
#' @export
#'
#' @examples
setGeneric("remove_sample", function(object, samples) standardGeneric("remove_sample"))

#' Title
#'
#' @param ProtData
#'
#' @return
#' @export
#'
#' @examples
setMethod("remove_sample",
          "ProtData",
          function(object, samples){
            dat <- object@data
            cond <- object@condition

            dat <- dat[,!(colnames(dat) %in% samples),drop = FALSE]
            cond <- cond[!rownames(cond) %in% samples, ,drop = FALSE]

            object@data <- dat
            object@condition <- cond
            return(object)
          }
)
