# graph.virome2.R
#'
#' Convert virome.df into an Virome igraph object - all nodes
#'
#' @param virome.df  data.frame, produced by the get.virome() function
#' @param no.stats   logical, do not calculate statistical tests [FALSE]
#' @param node1      string, column in virome.df to use as nodes ["run"]
#' @param node2      string, column in virome.df to use as edges ["sotu"]
#' @param con        pq-connection, use SerratusConnect()
#' @return graph.virome, a bipartite igraph object
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
                           no.stats   = FALSE,
                           con   = SerratusConnect()) {
  # Initialize Graph ================================================
  # Virome Graph objects by definition are run - sotu nodes
  node1 = 'run'
  node2 = 'sotu'
  
  # Create the edge List directly from 2 columns of virome.df
  edgeList <- virome.df[,  c(node1, node2)]

  # Create bipartite igraph object from Edge List
  g <- graph_from_data_frame(edgeList, directed = FALSE)
  
  # Color Nodes by node.type
  node1.match     <- unique( match(virome.df[ , node1 ], attributes(V(g))$names) ) # run
  node2.match     <- unique( match(virome.df[ , node2 ], attributes(V(g))$names) ) # sOTU
  
  # Define Vertex Types to make it bipartite
  V(g)$type               <- FALSE  # Run  nodes
  V(g)$type[node2.match]  <- TRUE   # sOTU nodes
  
  # Paint sOTU Nodes ================================================
    # Initialize sOTU / run data.frame for painting nodes
  sotu.df <- distinct( virome.df[ , c('sotu', 'nickname', 'gb_pid', 'gb_acc', 'tax_species', 'tax_family')])
  
  # set sOTU not mappted to GenBank to 0
  sotu.df$gb_pid[ sotu.df$gb_pid == -1 ] <- 0
  
  # Populate virome.df data into virGraph2
  node2.smatch <- match( attributes(V(g))$names ,  sotu.df$sotu )
    node2.smatch <- node2.smatch[ !is.na(node2.smatch)] 
  
  V(g)$sotu     <- "NA"
  V(g)$sotu[node2.match]        <- as.character( sotu.df$sotu[ node2.smatch ])
  
  V(g)$nickname <- "NA"
  V(g)$nickname[node2.match]    <- as.character( sotu.df$nickname[ node2.smatch ])
  
  V(g)$tax_speces <- "NA"
  V(g)$tax_species[node2.match] <- as.character(sotu.df$tax_species[ node2.smatch ])
  
  V(g)$tax_family <-  "NA"
  V(g)$tax_family[node2.match]  <- as.character(sotu.df$tax_family[ node2.smatch ])
  
  V(g)$gb_pid    <-  0
  V(g)$gb_pid[node2.match]     <- sotu.df$gb_pid[ node2.smatch ]
  
  V(g)$gb_acc    <-  "NA"
  V(g)$gb_acc[node2.match]     <- as.character(sotu.df$gb_acc[ node2.smatch ])
  
  # Paint Run Nodes =================================================
  run.df  <- distinct( virome.df[ , c('run', 'scientific_name', 'bio_sample', 'bio_project')])
  
  node1.rmatch <- match( attributes(V(g))$names,  run.df$run )
    node1.rmatch <- node1.rmatch[ !is.na(node1.rmatch)] 
  
  V(g)$run     <- "NA"
  V(g)$run[node1.match]        <- as.character( run.df$run[ node1.rmatch ])
  
  V(g)$scientific_name     <- "NA"
  V(g)$scientific_name[node1.match]        <- as.character( run.df$scientific_name[ node1.rmatch ])
  
  V(g)$bioproject     <- "NA"
  V(g)$bioproject[node1.match]        <- as.character( run.df$bio_project[ node1.rmatch ])
  
  # Paint Edges =====================================================
  # For visualization, it's useful to color edges by sOTU-node / run-node properties
  # sOTU properties
  E(g)$sotu            <- virome.df$sotu
  E(g)$node_coverage   <- virome.df$node_coverage # node == sotu here
  E(g)$tax_species     <- virome.df$tax_species
  E(g)$tax_family      <- virome.df$tax_family
  E(g)$gb_pid          <- virome.df$gb_pid
  
  # Run properties
  E(g)$run             <- virome.df$run
  E(g)$bioproject      <- virome.df$bio_project
  E(g)$scientific_name <- virome.df$scientific_name
  
  #Virome Enrichment Statistics ====================================
  if (no.stats){
    # skip statistics calculations
    V(g)$vrich     <-   -1
    V(g)$v.exact   <-   -1
    V(g)$v.or      <-    0
    V(g)$component <-    0
    V(g)$pr        <-    0
    V(g)$vrank     <-    0
    V(g)$lpa.label <-    NA
    
  }  else {
    # VRICH Calculation ---------------------------------
    # Get SRA-wide counts for each sotu-run
    sotu.gcount <-     tbl(con, "palm_virome_count") %>%
      dplyr::filter( sotu %in% sotu.df$sotu ) %>%
      select(sotu, runs) %>%  
      as.data.frame()
      dbDisconnect( con ) # Close Connection
      colnames(sotu.gcount) <- c("sotu", "n_total")
    
    # Get virome-wide counts for each sotu (same as degree)
    sotu.vcount <- data.frame( table(virome.df$sotu[ !duplicated(virome.df[c('run', 'sotu')]) ] ))
      colnames(sotu.vcount) <- c("sotu", "n_vir")
    
    # Merge dataframes, and calculate "vrich"
    # Percent sOTU in Virome vs SRA wide
    sotu.df <- merge( sotu.df, sotu.gcount, by = 'sotu')
    sotu.df <- merge( sotu.df, sotu.vcount, by = 'sotu')
    sotu.df$n_out <- ( sotu.df$n_total - sotu.df$n_vir )
    
    sotu.df$vrich <- sotu.df$n_vir / sotu.df$n_total
    
    # Assign
    node2.smatch <- match( attributes(V(g))$names ,  sotu.df$sotu )
    node2.smatch <- node2.smatch[ !is.na(node2.smatch)] 
    
    V(g)$vrich     <-   -1
    V(g)$vrich[node2.match]      <- sotu.df$vrich[ node2.smatch ]
    
    # VExact Calculation ------------------------------
    #  {Virome}   is the sum observation of all i sOTU in given virome
    #  {Serratus} is the sum observation of all sOTU across all Serratus
    
    #  N       : Count of all sOTU observations across {Virome}
    #  n_vir   : The observed ith sOTU count in {Virome} 
    #  n_out   : The observed ith sOTU count in {Serratus}, outside {V}
    #  n_total : The total count of ith sOTU observations
    
    #  M       : Count of all sOTU observations across {Serratus}
    #  m_vir   : The sum of all non-ith sOTU in {Virome}
    #  m_out   : The sum of all non-{V} sOTU in {Serratus}
    #  m_total : The total count of all non-ith sOTU observations
  
    
    # Take row-count from table meta-data (faster than counting rows)
    # SELECT reltuples::bigint AS m_total FROM pg_class WHERE relname = 'palm_virome';
    M_ <- 4076288
    N_ <- sum(sotu.df$n_vir)
    n.tests_ <- length(sotu.df$sotu)
    
    VExact <- function( sotu.n,
                        N = N_,
                        M = M_,
                        n.tests = n.tests_){
      # Return p-value of Fisher's Exact Test based on the above values
      # Operates on a single row
      n_vir   <- as.numeric( sotu.n[1] )
      n_total <- as.numeric( sotu.n[2] )
      n_out   <- n_total - n_vir
      
      m_vir <- N - n_vir
      m_out <- M - m_vir
      
      m_total <- M - n_total
      
      if ( FALSE ){
        # Set to True to print contigency table
        print(paste0(" Number of total tests: ",     n.tests ))
        print(paste0(" N : ",     N ))
        print(paste0(" M : ",     M ))
        print(paste0("    |  sOTU  |  !sOTU  ")              )
        print(paste0("{V} |  ", n_vir, " | ", m_vir )        )
        print(paste0("{S} |  ",  n_out, " | ", m_out )       )
        print(paste0("sum |  ",  n_total, " | ", m_total )   )
      }
        
      # fisher.test(rbind(c(1,9),c(11,3)), alternative="less")$p.value
      FT <- fisher.test( rbind(  c( n_vir, m_vir ),
                                 c( n_out, m_out )) ,
                         alternative = 'greater' )

      # Virome Exact Score w/ Bonferonni correction
      v.exact <- -log10( min(1 ,  FT$p.value * n.tests) )

      # IF Virome Exact is >100 or INF, set to 10
      if ( v.exact > 10 ){
        v.exact <- 10
      }
      
      # IF Virome Exact is == 0, set it to 0.1
      if ( v.exact == 0){
        v.exact <- 0.1
      }

      # IF odds ratio is "Infinite", set it to 10
      if ( is.infinite(FT$estimate) ){
        FT$estimate <- 10
      }

        # Return p-value, OR, and BJ corrected p-value as -log10
        return( c( FT$p.value,
                   v.exact,
                   FT$estimate)
                )
      
      # FAKE RETURN
      # TODO why is this error coming up
      # Error in fisher.test(rbind(c(n_vir, m_vir), c(n_out, m_out)), alternative = "greater") : 
      # all entries of 'x' must be nonnegative and finite
      # return( c( 0.05,
      #            1,
      #            1)
      # )
      
    }
    
    sotu.exact <- data.frame( t( apply(sotu.df[ ,  c("n_vir", "n_total")], VExact, MARGIN = 1 ) ) )
      colnames(sotu.exact) <- c("p.exact", "v.exact", "v.or")
  
    # Assign
    V(g)$v.exact   <-   -1
    V(g)$v.exact[node2.match]    <- sotu.exact$v.exact[ node2.smatch ]
    
    V(g)$v.or      <-   0
    V(g)$v.or[node2.match]       <- sotu.exact$v.or[ node2.smatch ]
    
    # Virome Rank ----------------------------------------
    # Calculate page_rank & virome_rank
    # Project node-type 2 (sOTU) into a monopartite network
    # weighted by number of observations
    mono.g <- bipartite_projection( g, multiplicity = TRUE, which = "true")
    
    # Calculate node2 page_rank
    V(mono.g)$pr <- page_rank( mono.g, weights = mono.g$weight )$vector
    V(g)$pr <- 0
    V(g)$pr[node2.match] <- V(mono.g)$pr
    
    # Assign
    # Calculate node2 virome_rank ( product of vrich and page_rank) 
    V(g)$vrank <- 0
    V(g)$vrank <- V(g)$pr * V(g)$vrich * V(g)$v.exact
    
    # Virome Community Detection ======================================
    # Components  -------------------------------------------
    comp.g <- components(g)
      Vc.label <- factor( comp.g$membership ) #original labels
    
      # Rename component-labels (membership) based on size
      # 1 = largest ... n = smallest
      rank.order <- order(comp.g$csize, decreasing = TRUE)
        levels(Vc.label) <- order(rank.order)
    
      # Assign labels to vertices
      V(g)$component <- as.character( Vc.label )
      
      # Label Propagation Algorithm -------------------------
      lpa.vrank <- function(g, labelname = 'scientific_name') {
  
        # Transfer Node Value (Vrank) to Edge (via sOTU)
        vir.nodes    <- data.frame( sotu  = V(g)$sotu,
                                    vrank = V(g)$vrank)
        vir.nodes  <- vir.nodes[ vir.nodes$sotu != "NA", ]
        vir.select <- match(levels( E(g)$sotu ), vir.nodes$sotu ) # ordered virus node assignment in g
  
        E(g)$vrank <- E(g)$sotu
          levels( E(g)$vrank ) <- vir.nodes$vrank[ vir.select ]
        E(g)$vrank <- as.numeric( as.character( E(g)$vrank ) )
  
        # Label Propagation Algorithm - Weighted by VRank
        lab.df <- data.frame( names  = V(g)$name,
                              component = V(g)$component,
                              label  = factor( vertex_attr( g, labelname ) ),
                              int.label  = -1 ,
                              fixed  = !V(g)$type )
        

        if (levels( lab.df$label )[1] == "NA" & length( levels( lab.df$label) ) == 2 ){
          # Edge Case: Convert labels to integers
          # NA is assigned integer value of 1, and other label is assigned 2
          # does not work with LPA, switch it
          lab.df$int.label[ lab.df$fixed ] <- 1
        } else {
          # Use Factor-numeric as integer labels
          lab.df$int.label[ lab.df$fixed ] <- as.numeric( lab.df$label[ lab.df$fixed ] )
        }
        
        lpa <- cluster_label_prop(g,
                                  weights = E(g)$vrank,
                                  initial = lab.df$int.label,
                                  fixed   = lab.df$fixed)
        
        lab.df$lpa       <- membership(lpa)
  
        # Re-name LPA clusters per label grouips
        label_groups   <- levels( lab.df$label )
        label_groups <- label_groups[ label_groups != "NA" ]
        label_groups <- lab.df[ match(label_groups, lab.df$label), c("label", "lpa")]
  
        lab.df$lpa.label <- label_groups$label[ match(lab.df$lpa, label_groups$lpa) ]
  
        V(g)$lpa.label   <- lab.df$lpa.label

        return( V(g)$lpa.label )
      }
    
      V(g)$lpa.label <- lpa.vrank(g)
  }
    
  return(g)
}
