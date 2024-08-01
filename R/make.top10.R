#' make.top10
#'
#' Aggregate a factor/vector and return only the top N entries
#' with all other entries renamed "Other"
#' Returned levels are rank-ordered by abundance for plotting
#'
#' @param invec   factor, unordered and repeated entries
#' @param top.n   numeric, how many entries to return before grouping to "Other" [10]
#' @param rename  character, the string to return instead of ["Other"]
#' 
#' @return factor, of the same length as invec, but with n+1 entries renamed "Other"
#' 
#' @keywords palmid Serratus palmdb sOTU
#' @examples
#' 
#' ## R Code Example
#' 
#' testvec <- c("A", "A", "A", "B", "C", "C", "C", "C")
#' 
#' makeTop10( testvec, top.n = 2 )
#' 
#' @export
# Retrieve nickname from an input of a palm_ids
makeTop10 <- function( invec , top.n = 10, rename = "Other"){
  # Convert to character
  invec <- as.character(invec)
  
  # Count entries, return top-10
  t10.entries <- table( invec )
  t10.entries <- t10.entries[ rev( order(t10.entries) ) ]
  t10.entries <- rownames(t10.entries)[1:top.n]
  
  # Refactor input vector to only top-10 entries
  # and levels ordered by top-10 order
  invec2 <- invec
  invec2[ !(invec2 %in% t10.entries) ] <- rename
  invec2 <- factor(invec2,
                   levels = c(t10.entries, rename))
  
  return(invec2)
  
}
