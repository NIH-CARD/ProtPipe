#' A protData contructor for SomaScan Data
#'
#' @param adat
#'
#' @return
#' @export
#'
#' @examples
create_protdata_from_soma <- function(adat, condition = NULL, filter = TRUE) {
  if(filter){
    dat <- soma_all_output(adat)
    dat <- Buffer_filter(dat)
  }else{
    dat <- soma_sample_out(adat)
  }

  return(create_protdata(dat, condition, method = "SomaScan"))
}


soma_sample_out=function(DT){
  annoaa <<- SomaDataIO::getAnalyteInfo(DT)
  anno <- SomaDataIO::getAnalyteInfo(DT)%>%
    dplyr::filter(Organism == "Human") %>%
    dplyr::filter(Type == "Protein")
  DT=as.data.frame(DT)
  DT_dat=DT%>%
    dplyr::filter(grepl("Sample", SampleType, ignore.case = TRUE))

  #check for duplicated SampleID
  duplicate_ids <- DT_dat$SampleId[duplicated(DT_dat$SampleId) | duplicated(DT_dat$SampleId, fromLast = TRUE)]
  if(length(duplicate_ids>0)){
    cat(paste0("removing duplicates: ", paste(duplicate_ids, collapse = ", ")))
    DT_dat <- DT_dat[!DT_dat$SampleId %in% duplicate_ids, ]
  }


  rownames(DT_dat)=DT_dat$SampleId
  DT_dat=DT_dat%>%
    dplyr::select(matches("seq\\.", ignore.case = TRUE))%>%
    t()
  DT_dat=merge(anno[,grep('AptName|UniProt|EntrezGeneSymbol|TargetFullName',colnames(anno))], DT_dat,by.x='AptName',by.y=0,all.x=T)
  DT_dat=DT_dat %>%
    dplyr::filter(UniProt != "")%>%
    dplyr::filter(EntrezGeneSymbol != "") %>%
    dplyr::rename(Protein_Group= UniProt)%>%
    dplyr::rename(Genes= EntrezGeneSymbol)
  return(DT_dat)
}

##format data
soma_all_output=function(DT){
  anno=SomaDataIO::getAnalyteInfo(DT)
  DT=data.frame(DT)
  DT_dat=data.frame(DT)%>%
    dplyr::filter(grepl("Sample", SampleType, ignore.case = TRUE))

  #check for duplicated SampleID
  duplicate_ids <- DT_dat$SampleId[duplicated(DT_dat$SampleId) | duplicated(DT_dat$SampleId, fromLast = TRUE)]
  if(length(duplicate_ids>0)){
    cat(paste0("removing duplicates: ", paste(duplicate_ids, collapse = ", ")))
    DT_dat <- DT_dat[!DT_dat$SampleId %in% duplicate_ids, ]
  }

  rownames(DT_dat)=DT_dat$SampleId
  DT_dat = DT_dat[, grep('seq\\.', colnames(DT_dat))] %>%
    t() %>%
    as.data.frame()

  Buffer_mean <- DT %>%
    dplyr::filter(SampleType == "Buffer") %>%
    dplyr::select(matches("seq\\.", ignore.case = TRUE)) %>%
    dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    t() %>%
    as.data.frame()%>%
    {colnames(.) <- 'Buffer'; .}

  Calibrator_mean <- DT %>%
    dplyr::filter(SampleType == "Calibrator") %>%
    dplyr::select(matches("seq\\.", ignore.case = TRUE)) %>%
    dplyr::summarise(across(everything(), ~ mean(.x, na.rm = TRUE))) %>%
    t() %>%
    as.data.frame()%>%
    {colnames(.) <- 'Calibrator'; .}

  DT_combined <- cbind(Buffer_mean, Calibrator_mean, DT_dat)
  DT_out=merge(anno[,grep('AptName|UniProt|EntrezGeneSymbol|TargetFullName|Organism|Type',colnames(anno))], DT_combined,by.x='AptName',by.y=0)
  DT_out=DT_out%>%
    dplyr::rename(Protein_Group= UniProt)%>%
    dplyr::rename(Genes= EntrezGeneSymbol)
  return(DT_out)
}

Buffer_filter=function(DT){
  DT=as.data.frame(DT)
  DT_filter <<- DT %>%
    dplyr::mutate(across(
      .cols = -c(Protein_Group, Genes, Buffer,Calibrator),  # Exclude PG_group, genes, and Buffer
      .fns = ~ ifelse(. < Buffer, NA, .)  # Apply the condition
    ))
  return(DT_filter)
}
