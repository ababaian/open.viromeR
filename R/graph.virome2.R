# graph.virome2.R
#'
#' Convert virome.df into an igraph object - all node
#'
#' @param virome.df  data.frame, produced by the get.virome() function
#' @param node1      string, column in virome.df to use as nodes ["bio_project"]
#' @param node2      string, column in virome.df to use as edges ["tax_family"]
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
graph.virome2  <- function(virome.df  = NA,
                           node1 = 'run',
                           node2 = 'sotu',
                           node1.col = 'bio_project',
                           node2.col = 'tax_family',
                           edge.col = 'node_coverage') {
  
  # Create the edge List directly from 2 columns of virome.df
  edgeList <- virome.df[,  c(node1, node2)]

  # Create igraph object from Edge List
  g <- graph_from_data_frame(edgeList, directed = FALSE)
  
  # Paint edges direct from virome.df
  E(g)$sotu          <- virome.df$sotu
  E(g)$run           <- virome.df$run
  E(g)$node_coverage <- virome.df$node_coverage
  E(g)$gb_pid        <- virome.df$gb_pid
  E(g)$bioproject    <- virome.df$bio_project
  E(g)$scientific_name <- virome.df$scientific_name
  
  # Color Nodes by node.col
  node1.match     <- unique( match(virome.df[ , node1 ], attributes(V(g))$names) )
  node2.match     <- unique( match(virome.df[ , node2 ], attributes(V(g))$names) )
  
  V(g)$layer              <- node2
  V(g)$layer[node1.match] <- node1
  
  # Paint sotu nodes
  if (node2 == 'sotu'){
    sotu.df <- distinct( virome.df[ , c('sotu', 'nickname', 'gb_pid', 'gb_acc', 'tax_species', 'tax_family')])
    
    # Get SRA-wide counts for each sotu
    sotu.gcount <-     tbl(con, "palm_virome_count") %>%
      dplyr::filter( sotu %in% sotu.df$sotu ) %>%
      select(sotu, runs) %>%  
      as.data.frame()
    # Get virome-wide counts for each sotu
    sotu.vcount <- data.frame( table(virome.df$sotu) )
    
    # Merge, and calculate "vrich"
    # Percent sOTU in Virome vs SRA wide
    sotu.df <- merge( sotu.df, sotu.gcount, by = 'sotu')
    sotu.df <- merge( sotu.df, sotu.vcount, by.x = 'sotu', by.y = 'Var1')
    sotu.df$vrich <- round( 100 * sotu.df$Freq / sotu.df$runs, 2)
    # set unmapped to 0
    sotu.df$gb_pid[ sotu.df$gb_pid == -1 ] <- 0
    
    
    # Populate virome.df data into virGraph2
    node2.smatch <- match( attributes(V(g))$names ,  sotu.df$sotu )
    node2.smatch <- node2.smatch[ !is.na(node2.smatch)] 
    
    V(g)$sotu     <- "NA"
    V(g)$sotu[node2.match]      <- as.character( sotu.df$sotu[ node2.smatch ])
    
    V(g)$nickname <- "NA"
    V(g)$nickname[node2.match]  <- as.character( sotu.df$nickname[ node2.smatch ])
    
    V(g)$tax_speces <- "NA"
    V(g)$tax_species[node2.match] <- as.character(sotu.df$tax_species[ node2.smatch ])
    
    V(g)$tax_family <-  "NA"
    V(g)$tax_family[node2.match] <- as.character(sotu.df$tax_family[ node2.smatch ])
    
    V(g)$gb_pid    <-  -1
    V(g)$gb_pid[node2.match]     <- sotu.df$gb_pid[ node2.smatch ]
    
    V(g)$gb_acc    <-  "NA"
    V(g)$gb_acc[node2.match] <- as.character(sotu.df$gb_acc[ node2.smatch ])
    
    V(g)$vrich     <-   101
    V(g)$vrich[node2.match] <- sotu.df$vrich[ node2.smatch ]
    
  }
  
  # Calculate Components (communities) of the graph
  V(g)$component <- components(g)$membership
  
  return(g)
}
