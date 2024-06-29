#' get.sraProj
#'
#' Retrieve the "BioProject" field for a set of SRA 'run_id'
#' and return all SRA run_id for each unique BioProject
#' <UNDER DEVELOPMENT>
#' 
#' @param run_ids character, SRA 'run_id'
#' @param con     pq-connection, use SerratusConnect()
#' @param exclude.input.runs boolean, exclude runs input 'run_ids'
#' @return data.frame, run_id, biosample character vectors
#' @keywords palmid Serratus BioProject SRA
#' @examples
#' \donttest{
#' # SRA Library of interest
#' con <- SerratusConnect()
#' library.bioProject   <- get.sraProj( 'SRR9968562' , con)
#' }
#' @import RPostgreSQL
#' @import dplyr ggplot2
#' @export
get.sraProj <- function(run_ids, con, exclude.input.runs = FALSE) {
  # Bind Local Variables
  run <- bio_sample <- biosample_id <-bioproject <- NULL
  
  # Forward Search to retrieve BioProjects
  # get biosample field for run_id
  sra.bio <- tbl(con, "srarun") %>%
    dplyr::filter(run %in% run_ids) %>%
    select(run, bio_sample, bio_project) %>%
    as.data.frame()
    colnames(sra.bio) <- c("run_id", "biosample_id", "bioproject")
    # must be unique
    sra.bio <- sra.bio[ !duplicated(sra.bio$run_id), ]
    bp_ids  <- sra.bio$bioproject [ !duplicated(sra.bio$bioproject) ]
    
  # Reverse Search to retrieve all runs in BioProject
  sra.bio2 <- tbl(con, "srarun") %>%
      dplyr::filter(bio_project %in% bp_ids) %>%
      select(run, bio_sample, bio_project) %>%
      as.data.frame()
    colnames(sra.bio2) <- c("run_id", "biosample_id", "bioproject")
    # must be unique
    sra.bio2 <- sra.bio2[ !duplicated(sra.bio2$run_id), ]
    
  # Deplete runs in the input set of 'run_ids'
  if (exclude.input.runs) {
    sra.bio2 <- sra.bio2[ !(sra.bio2$run_id %in% sra.bio$run_id), ]
  }
    
    return(sra.bio2)
  # 
}
