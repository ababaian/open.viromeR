# graph.palm.R
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
graph.palm     <- function(virome.df = NA,
                           pid.threshold = 30,
                           con = SerratusConnect()) {
  # Reduce sOTU list to unique entries
  # sotu.vec <- unique(virome.df$sotu)
  sotu.vec <- as.character( unique(virome.df$sotu) )
  
  # Get palmVirome based on sra.vec
  # add self-identity (100%)
  sotu.edge <- tbl(con, "palm_graph") %>%
    dplyr::filter( palm_id1 %in% sotu.vec ) %>%
    as.data.frame()
  
  # sotu.edge <- rbind(sotu.edge,
  #                    data.frame(row_index = 1,
  #                               palm_id1 = sotu.vec,
  #                               palm_id2 = sotu.vec,
  #                               pident = 100) )
  sotu.edge <- sotu.edge[ sotu.edge$palm_id2 %in% sotu.vec,  ]
  sotu.edge <- sotu.edge[ sotu.edge$pident >= pid.threshold, ]
  
  # Paint vertices/edges from virome.df
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
    sotu.df$vrich <- round( sotu.df$Freq / sotu.df$runs, 2)
    sotu.df <- sotu.df[,  c('sotu', 'nickname', 'gb_pid', 'gb_acc', 'tax_species', 'tax_family', 'vrich')]
  
  # Create igraph object from Edge List
  g <- graph_from_data_frame(sotu.edge[, c("palm_id1", "palm_id2")],
                             directed = FALSE,
                             vertices = sotu.df)
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
  
  return(g)
}