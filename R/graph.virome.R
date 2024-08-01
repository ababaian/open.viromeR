# graph.virome.R
#'
#' Convert virome.df into an igraph object
#'
#' @param virome.df  data.frame, produced by the get.virome() function
#' @param nodes      string, column in virome.df to use as nodes ["bio_project"]
#' @param edges      string, column in virome.df to use as edges ["tax_family"]
#' @param edge.threshold numeric, minium number of edges required to report [1]
#' @param node.col   string, column in virome.df to use for node colors ["bio_project"]
#' @param edge.col   string, column in virome.df to use for edge colors ["sotu"]
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
graph.virome   <- function(virome.df  = NA,
                           nodes = 'bio_project',
                           edges = 'tax_family',
                           edge.threshold = 1,
                           node.col = 'bio_project',
                           edge.col = 'sotu') {
  
  # Create the nodes and edge Table
  virTable <- table( virome.df[,  c(edges, nodes)] )
  
  # Create a network
  row2edgeList <- function(input.row, n.threshold = edge.threshold){
    # For each edge-type, retrieve edgeList of between nodes
    edgeList <- names(input.row[ which(input.row >= n.threshold) ])
    
    if ( length(edgeList) <= 1){
      # No edges of this type are found
      return(NULL)
    } else {
      # Generate edgeList for input row, and recipricol
      # by taking all combinations of nodes with an edge-type
      edgeList <- matrix( c(combn(edgeList, 2, t)), byrow = T, ncol = 2)
      # Calculate recipricol (for directed graphs)
      #edgeList <- rbind( edgeList[, c(1,2)], edgeList[, c(2,1)])
      return(edgeList)
    }
  }
  
  # A List of Edgelists per type of edges
  edgeList.g   <- apply( virTable, 1, row2edgeList)
  
  # Initialize Edgelist Matrix
  n = 1
  while ( is.null(edgeList.g[[n]]) ){
    n <- n+1
  }
  edgeList.mat <- matrix( unlist( edgeList.g[n] ), ncol = 2, byrow = F)
  rownames(edgeList.mat) <- rep( names(edgeList.g)[n], length(edgeList.mat)/2)
  
  # Populate EdgeList Matrix
  for (n in (n+1):length(edgeList.g)){
    if ( !is.null(edgeList.g[[n]]) ){
      edgeList.tmp <- matrix( unlist( edgeList.g[n] ), ncol = 2, byrow = F)
      rownames(edgeList.tmp) <- rep( names(edgeList.g)[n], length(edgeList.tmp)/2)
    
      edgeList.mat <- rbind(edgeList.mat,
                            edgeList.tmp)
    }
  }
  
  # Create weighted adjacency matrix sum of all edge-types
  #adjMat <- as_adjacency_matrix( graph_from_data_frame(edgeList.mat, directed = FALSE) )
  #
  # Create igraph object from Adj Matrix
  #g      <- graph_from_adjacency_matrix( adjMat, weighted = 'Weights')
  #  E(g)$arrow.mode <- "-"
  #  E(g)$Weights <- E(g)$Weights * (10/max(E(g)$Weights)) # scale to max 10 weight
  
    
  # Create igraph object from Edge List
  g <- graph_from_data_frame(edgeList.mat, directed = FALSE)
    
  # Color Nodes by node.col
  node.match      <- match(attributes(V(g))$names,  virome.df[ , nodes ])
  V(g)$node.col   <- virome.df[ node.match, node.col ]
  V(g)$vcol       <- as.character(plyr::mapvalues(V(g)$node.col,
                                                  from = levels(V(g)$node.col),
                                                  to = rainbow( length( levels(V(g)$node.col)))))
  
  # Color Edges by edge.col
  #edge.match      <- match(attributes(E(g))$vnames,  virome.df[ , edges ])
  E(g)$edge.col   <- rownames(edgeList.mat)
  E(g)$ecol       <- as.character(plyr::mapvalues(E(g)$edge.col,
                                                  from = levels(E(g)$edge.col),
                                                  to = rainbow( length( levels(E(g)$edge.col)))))
  
  # Calculate Components (communities) of the graph
  V(g)$component <- components(g)$membership
  
  return(g)
}

