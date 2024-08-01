# graph.palm.R
#'
#' Convert virome.df into a palmprint-palmprint igraph object
#'
#' @param virome.df   data.frame, produced by the get.virome() function
#' @param sotu.vec    character, vector of sOTU id to retrieve
#' @param pid.threshold numeric, min percent identity to include an edge [40]
#' @param con        pq-connection, use SerratusConnect()
#' @return graph.virome, an igraph object
#' @keywords palmid sql sra biosample bioproject Serratus igraph
#' @examples
#'
#' \donttest{
#' graph.virome(virome.df, nodes = "run", edges = "sotu", node.col = "bio_project")
#' }
#'
#' @import dplyr igraph
#' @export
#'
graph.palm     <- function(virome.df = NA,
                           sotu.vec  = NA,
                           pid.threshold = 30,
                           expanded.graph = FALSE,
                           all.edges      = FALSE,
                           con = SerratusConnect()) {
  
  if ( !is.na( virome.df )[1] ){
    # virome.df provided
    # Paint vertices/edges from virome.df
    sotu.df <- distinct( virome.df[ , c('sotu', 'nickname', 'gb_pid', 'gb_acc', 'tax_species', 'tax_family')])
    # Reduce sOTU list to unique entries
    sotu.vec <- as.character( sotu.df$sotu )
  } else if ( !is.na( sotu.vec ) ) {
    sotu.vec <- as.character( sotu.vec )
    sotu.df  <- data.frame( sotu = sotu.vec,
                            nickname = NA,
                            gb_pid   = NA,
                            gb_acc   = NA,
                            tax_species = NA,
                            tax_family  = NA)
  } else {
    stop("virome.df or sotu.vec must be provided as input.")
  }
    
  # Get palmVirome based on sra.vec
  # add self-identity (100%)
  sotu.edge <- tbl(con, "palm_graph") %>%
    dplyr::filter( palm_id1 %in% sotu.vec ) %>%
    as.data.frame()
    dbDisconnect( con ) # Close Connection
  
  if (expanded.graph){
    # Return all linked sOTU
    sotu.edge <- sotu.edge[ sotu.edge$pident >= pid.threshold, ]
    
    # Create igraph object from Edge List
    g <- graph_from_data_frame(sotu.edge[, c("palm_id1", "palm_id2")],
                               directed = FALSE)
    E(g)$pid <- sotu.edge$pident
    
    if (!all.edges){
      # Prune edges to those with Degree >= 2
      # (removes single-link viruses)
      g <- subgraph( g, vids = which(degree(g) >= 2) )
    }
    
  } else {
    # Return sOTU in input set only
    sotu.edge <- sotu.edge[ sotu.edge$palm_id2 %in% sotu.vec,  ]
    sotu.edge <- sotu.edge[ sotu.edge$pident >= pid.threshold, ]
    
    # Create igraph object from Edge List
    g <- graph_from_data_frame(sotu.edge[, c("palm_id1", "palm_id2")],
                               directed = FALSE,
                               vertices = sotu.df)
    E(g)$pid <- sotu.edge$pident
  }
  
  # Calculate Components (communities) of the graph
  comp.g <- components(g)
  Vc.label <- factor( comp.g$membership ) #original labels
  
  # Rename component-labels (membership) based on size
  # 1 = largest ... n = smallest
  rank.order <- order(comp.g$csize, decreasing = TRUE)
  levels(Vc.label) <- order(rank.order)
  
  # Assign labels to vertices
  V(g)$component <- as.character( Vc.label )
  
  return(g)
}
