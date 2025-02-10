#' SILAC Class
#'
#' An S4 class that holds SILAC proteomics data and provides methods for processing.
#'
#' @export
setClass("SilacData",
         contains = "ProtData",
         slots = c(
           channels = "list"
          )
         )


create_silac <- function(dat, intensity_cols, channel_names = c("Channel1", "Channel2"), condition = NULL, method = "SILAC") {

  #create ProtData object
  pdata <- create_protdata(dat = dat, intensity_cols = intensity_cols, condition = condition, method = method)

  # Split the data into the individual channels
  channels <- list()
  for (channel in channel_names) {
    # Get columns corresponding to the current channel
    channel_cols <- grep(channel, colnames(dat), value = TRUE)
    if (length(channel_cols) == 0) {
      warning(paste("No columns found for channel:", channel))
      channels[[channel]] <- NULL
    } else {
      channel_data <- dat[, channel_cols, drop = FALSE] # Subset channel columns
      colnames(channel_data) <- trim_colnames(channel_data) # Trim column names
      channels[[channel]] <- channel_data
    }
  }

  #return final object
  new("SilacData",
      data = pdata@data,              # Inherited from ProtData
      condition = pdata@condition,   # Inherited from ProtData
      prot_meta = pdata@prot_meta,   # Inherited from ProtData
      method = pdata@method,         # Inherited from ProtData
      channels = channels                  # Additional slot
  )
}

