# graph.palm.R
#'
#' Convert virome.df into a palmprint-palmprint igraph object
#'
#' @param sotu.vec       character, vector of sOTU id to retrieve
#' @param pid.threshold  numeric, min percent identity to include an edge [30]
#' @param expanded.graph logical, return non-input sOTU joined to >=2 input sOTU
#' @param all.edges      logical, return all sOTU with alignment to input sOTU
#' @param get.metadata   logical, return sOTU meta-data for all returned entries
#' @param con        pq-connection, use SerratusConnect()
#' @return graph.virome, an igraph object
#' @keywords palmid sql sra biosample bioproject Serratus igraph
#' @examples
#'
#' \donttest{
#' graph.palm(virome.df$sotu)
#' }
#'
#' @import dplyr igraph
#' @export
#'
graph.palm     <- function(sotu.vec       = NA,
                           pid.threshold  = 30,
                           expanded.graph = FALSE,
                           all.edges      = FALSE,
                           get.metadata   = TRUE,
                           con = SerratusConnect()) {
  
  # SOTU.VEC INPUT --------------------
  if ( !is.na( sotu.vec )[1] ) {
    sotu.vec <- unique( as.character( sotu.vec ) )
  } else {
    stop("sotu.vec must be provided as input.")
  }
  
  # Get palmVirome based on sra.vec
  # add self-identity (100%)
  sotu.edge <- tbl(con, "palm_graph") %>%
    dplyr::filter( palm_id1 %in% sotu.vec ) %>%
    as.data.frame()
  
  # EXPANDED PALMNET ------------------
  if (!expanded.graph){
    # Return sOTU in input set only
    sotu.edge <- sotu.edge[ sotu.edge$palm_id2 %in% sotu.vec,  ]
  }

  # MERGE PAIRWISE ALIGNMENTS -----------
  mergePairwiseAlignment <- function(sotu.edge_ = sotu.edge){
    # For each sOTU_X and sOTU_Y, both X_Y and Y_X alignments
    # are reported, average those pid and report single
    # alignment
    sotu.edge_ = sotu.edge
    id1LT <- (sotu.edge$palm_id1 < sotu.edge$palm_id2)
    
    sotu.edge2 <- sotu.edge_
    sotu.edge2$id1LT <- id1LT
    
    # Transpose palm_id1 and palm_id2 when 1 < 2
    sotu.edge2$palm_id1[ id1LT ]  <- sotu.edge_$palm_id2[ id1LT ]
    sotu.edge2$palm_id2[ id1LT ]  <- sotu.edge_$palm_id1[ id1LT ]
    
    sotu.edge2 <- aggregate(pident ~ palm_id1 + palm_id2, sotu.edge2, mean)
    
    return( sotu.edge2 )
  }
  sotu.edge <- mergePairwiseAlignment()
  
  # PID FILTER EDGES ---------------------
  sotu.edge <- sotu.edge[ sotu.edge$pident >= pid.threshold, ]
  
  # GET METADATA ------------------------
  sotu.nodes  <- unique( c( sotu.edge$palm_id1, sotu.edge$palm_id2 ))
  
  if (get.metadata){
    # Retrieve Meta-data
    sotu.df <- tbl(con, "palm_tax") %>%
      dplyr::filter( sotu %in% sotu.nodes ) %>%
      dplyr::filter( centroid == TRUE ) %>%
      select(sotu, nickname, percent_identity, gb_acc, tax_species, tax_family) %>%  
      as.data.frame()
      colnames(sotu.df) <- c("sotu", "nickname", "gb_pid", "gb_acc", "tax_species", "tax_family")
  } else {
    # NA for meta-data
    sotu.df <- data.frame( sotu = sotu.nodes,
                            nickname = NA,
                            gb_pid   = NA,
                            gb_acc   = NA,
                            tax_species = NA,
                            tax_family  = NA)
  }
  
  
  # Create igraph object from Edge List
  g <- graph_from_data_frame(sotu.edge[, c("palm_id1", "palm_id2")],
                             directed = FALSE,
                             vertices = sotu.df)
  E(g)$pid <- sotu.edge$pident
  
  # Graph Stats =========================================
  # Calculate Components (communities) of the graph
  comp.g <- components(g)
  Vc.label <- factor( comp.g$membership ) #original labels
  
  # Rename component-labels (membership) based on size
  # 1 = largest ... n = smallest
  rank.order <- order(comp.g$csize, decreasing = TRUE)
  levels(Vc.label) <- order(rank.order)
  
  # Assign labels to vertices
  V(g)$component <- as.character( Vc.label )
  
  dbDisconnect( con ) # Close Connection
  return(g)
}
