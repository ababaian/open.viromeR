# get.negativeVirome
#'
#' A wrapper of get* functions to create an empty palm.virome data.frame
#' which includes all runs analyzed in Serratus, not simply those with viruses
#' detected. Used to calculate how manmy "Virus Negative" datasets are in
#' a search.
#'
#' @param run.vec    character vector, sra run accession list
#' @param org.search string, search term to query "scientific_name". Uses SQL "like"
#' @param con       pq-connection, use SerratusConnect()
#' @return palm.virome  data.frame
#' @keywords palmid sql sra geo biosample bioproject timeline Serratus control
#' @examples
#'
#' con <- SerratusConnect()
#'
#' \donttest{
#' get.negativeVirome(org.search = 'Drosophilia%', con )
#' }
#'
#' @import dplyr ggplot2
#' @export
#'
get.negativeVirome <- function(run.vec    = NA,
                               org.search = NA,
                               con = SerratusConnect() ) {
  
  if ( !is.na(run.vec)[1] ){
    # Get palmVirome based on sra.vec 
    virome.df <- tbl(con, "srarun") %>%
      dplyr::select( run, scientific_name, bio_sample, bio_project) %>%
      dplyr::filter( run %in% run.vec ) %>%
      as.data.frame()
  } else if ( !is.na(org.search)[1] ) {
    # Get palmVirome based on SQL "like" search of org.search
    # agains the "scientific_name" column
    virome.df <- tbl(con, "srarun") %>%
      dplyr::select( run, scientific_name, bio_sample, bio_project ) %>%
      dplyr::filter( scientific_name %like% org.search ) %>%
      as.data.frame()
  } else {
    stop(paste0("One of 'sra.vec' or 'org.search' must be provided"))
  }
  
  # Check that return is non-empty
  if ( nrow(virome.df) == 0 ){
    stop(paste0("Error: no runs were returned."))
  } else {
    # Set each column class explicitly
    virome.df$run             <- factor(virome.df$run)
    virome.df$scientific_name <- factor(virome.df$scientific_name)
    virome.df$bio_sample      <- factor(virome.df$bio_sample)
    virome.df$bio_project     <- factor(virome.df$bio_project)
    virome.df$palm_id         <- NA 
    virome.df$sotu            <- NA 
    virome.df$nickname        <- NA
    virome.df$gb_acc          <- NA 
    virome.df$gb_pid          <- NA 
    virome.df$gb_eval         <- NA 
    virome.df$tax_species     <- NA
    virome.df$tax_family      <- NA
    virome.df$node            <- NA
    virome.df$node_coverage   <- NA 
    virome.df$node_pid        <- NA
    virome.df$node_eval       <- NA
    virome.df$node_qc         <- NA
    virome.df$node_seq        <- NA
  }
  
  # Add time (release date) to virome.df
  # virome.df$date <- get.sraDate(virome.df$run_id, con, TRUE)
  
  # Add geo-data if available to virome.df
  # virome.geo.tmp <- get.sraGeo( run_ids = NULL,
  #                               biosample_ids = virome.df$bio_sample,
  #                               con = con, ordinal =  TRUE)
  # 
  # if (!all(virome.geo.tmp$biosample_id == virome.df$biosample_id)) {
  #   stop("Error in geo lookup.")
  # } else {
  #   virome.df$lng <- virome.geo.tmp$lng
  #   virome.df$lat <- virome.geo.tmp$lat
  # }
  
  return(virome.df)
}

