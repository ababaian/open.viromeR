#' get.taxRunlist
#'
#' Given an input taxonomic term, retrive all SRA RunInfo runs matching that
#' term
#'
#' @param genus   character, exact character sequence
#' @param con     pq-connection, use SerratusConnect()
#' @return        data.frame, run_ids
#' @keywords palmid Serratus sra taxonomy
#' @examples
#' \donttest{
#' # SRA Library of interest
#' con <- SerratusConnect()
#' homo.list   <- get.taxRunlist( genus = 'Homo', con)
#' }
#' @import RPostgreSQL
#' @import dplyr ggplot2
#' @export
get.taxRunlist <- function(genus = NULL, con = SerratusConnect() ) {
  # Bind Local Variables
  tax.runs <- NULL

  # get run_id for input genus field for run_id
  tax.runs <- tbl(con, "sra_tax") %>%
    dplyr::filter(tax_genus == genus) %>%
    select(run) %>%
    as.data.frame()
    colnames(tax.runs) <- c("run")
    # must be unique
    tax.runs <- tax.runs[ !duplicated(tax.runs$run), ]
    
    DBI::dbDisconnect(con)
    return( tax.runs )
    
}
