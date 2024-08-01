#' get.sraMeta
#'
#' Retrieve the "BioSample", "BioProject", or "Both" field for a set of 
#' input SRA 'run_id'
#'
#' @param run_ids character, SRA 'run_id'
#' @param con     pq-connection, use SerratusConnect()
#' @param ordinal boolean, return 'run_ids' ordered vector [FALSE]
#' @return data.frame, SRA Run Metadata
#' @keywords palmid Serratus sra biosample bioproject
#' @examples
#' \donttest{
#' # SRA Library of interest
#' con <- SerratusConnect()
#' sra.df   <- get.sraMeta( 'SRR9968562' , con)
#' }
#' @import RPostgreSQL
#' @import dplyr ggplot2
#' @export
get.sraMeta <- function(run_ids, con, ordinal = FALSE) {
  # Bind Local Variables
  run <- release_date <- spots <- bio_project <- scientific_name <- library_strategy <- bases <- NULL
  
  # get sra metadata
  sra.meta <- tbl(con, "srarun") %>%
    dplyr::filter(run %in% run_ids) %>%
    select(run, release_date, spots, bases, bio_project, scientific_name, library_strategy ) %>%
    as.data.frame()
  colnames(sra.meta) <- c("run_id", "date", "spots", "bases", "bioproject", "scientific_name", "library_strategy")
  # must be unique
  sra.meta <- sra.meta[ !duplicated(sra.meta$run_id), ]
  
  if (ordinal){
    # Ordinal return
    # Left join on palm_ids to make a unique vector
    ord.sra <- data.frame( run_id = run_ids )
    ord.sra <- merge(ord.sra, sra.meta, all.x = TRUE)
    ord.sra <- ord.sra[ match(run_ids, ord.sra$run_id), ]
    
    return( ord.sra )
   } else {
    # Non-Ordinal, return unique entries
     return( sra.meta )
   }
}
