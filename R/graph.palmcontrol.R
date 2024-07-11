# graph.palmControl.R
#'
#' Convert virome.df into a palmprint-palmprint igraph object
#'
#' @param virome.df  data.frame, produced by the get.virome() function
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
graph.palmControl <- function(virome.df = NA,
                              pid.threshold = 30,
                              samples = 1,
                              con = SerratusConnect()) {
  # Extract list of all current sOTU
  sotu.all <- tbl(con, "palmdb2") %>%
    dplyr::filter( centroid == TRUE ) %>%
    select(sotu) %>%
    as.data.frame()
  
  # How many sOTU are in the input virome
  sotu.n <- length(unique(virome.df$sotu))
  
  for (i in 1:samples){
    # Sample N number of random sOTU from palmDB
    iset <- sample( sotu.all$sotu, sotu.n, replace = FALSE)
    
    # Get palmVirome based on sra.vec
    # add self-identity (100%)
    sotu.edge <- tbl(con, "palm_graph") %>%
      dplyr::filter( palm_id1 %in% iset ) %>%
      as.data.frame()
    DBI::dbDisconnect(con)
    
    sotu.edge <- sotu.edge[ sotu.edge$palm_id2 %in% iset,  ]
    sotu.edge <- sotu.edge[ sotu.edge$pident >= pid.threshold, ]
    
    # Create igraph object from Edge List
    g <- graph_from_data_frame(sotu.edge[, c("palm_id1", "palm_id2")],
                               directed = FALSE,
                               vertices = iset)
    E(g)$pid <- sotu.edge$pident
    
    # Calculate Components (communities) of the graph
    comp.g <- components(g)
    Vc.label <- factor( comp.g$membership ) #original labels
    
    # Rename component-labels (membership) based on size
    # 1 = largest ... n = smallest
    rank.order <- order(comp.g$csize, decreasing = TRUE)
    levels(Vc.label) <- order(rank.order)
    
    # Assign labels to vertices
    V(g)$component <- as.character( Vc.label )
    
  }
  return(g)
}
