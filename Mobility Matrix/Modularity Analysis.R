library(igraph)
library(ggplot2)

edge1=read.csv(file="class1to2_econ_node.csv")
edge1 = edge1[,-1]
head(edge1) 

edge2=read.csv(file="class2to3_econ_node.csv")
edge2 = edge2[,-1]
head(edge1) 

sem_edge1=read.csv(file="semester1to2_econ_node.csv")
sem_edge1 = sem_edge1[,-1]

sem_edge_di1=read.csv(file="semester1to2_econ_node_di.csv")
sem_edge_di1 = sem_edge_di1[,-1]


attr1=read.csv(file="class1to2_econ_attr.csv")
attr2=read.csv(file="class2to3_econ_attr.csv")
sem_attr1 = read.csv(file="semester1to2_econ_attr.csv")
sem_attr1 = sem_attr1[,-1]


attr1 = attr1[,-1]
attr2 = attr2[,-1]

head(attr1) 


edge_graph=graph_from_data_frame(d=edge1, directed=T,
                                    vertices=attr1) 
edge_graph 

edge_graph2=graph_from_data_frame(d=edge2, directed=T,
                                 vertices=attr2) 

edge_graph_sem=graph_from_data_frame(d=sem_edge1, directed=T,
                                 vertices=sem_attr1) 
edge_graph_sem=graph_from_data_frame(d=sem_edge_di1, directed=T,
                                     vertices=sem_attr1) 

plot(edge_graph, vertex.label=NA, edge.arrow.size=.5, 
     edge.arrow.width=.5, edge.color="light gray", 
     edge.curved=.2, vertex.frame.color= NA, margin=0)

plot(edge_graph2, vertex.label=NA, edge.arrow.size=.5, 
     edge.arrow.width=.5, edge.color="light gray", 
     edge.curved=.2, vertex.frame.color= NA, margin=0)

plot(edge_graph_sem, vertex.label=NA, edge.arrow.size=.5, 
     edge.arrow.width=.5, edge.color="light gray", 
     edge.curved=.2, vertex.frame.color= NA, margin=0)

edge_density(edge_graph)
edge_density(edge_graph2)
edge_density(edge_graph_sem)


cluster1=cluster_walktrap(graph=edge_graph, steps=4, membership=T)
cluster2=cluster_walktrap(graph=edge_graph2, steps=4, membership=T)
cluster1_sem=cluster_walktrap(graph=edge_graph_sem, steps=4, membership=T)
cluster1_sem_2=cluster_walktrap(graph=edge_graph_sem, steps=3, membership=T)


cluster1_m=modularity(cluster1)
cluster2_m=modularity(cluster2)

cluster1_sem_m=modularity(cluster1_sem)
cluster1_sem_m2=modularity(cluster1_sem_2)

cluster2_m
cluster1_sem_m

cluster1_c = membership(cluster1)
cluster2_c = membership(cluster2)
cluster2_sem_c = membership(cluster1_sem)
cluster2_sem_c2 = membership(cluster1_sem_2)


layout=layout.fruchterman.reingold(edge_graph) 
plot(edge_graph, layout=layout, #note the use of layout
     vertex.color=cluster1_c, edge.arrow.size=.5,
     edge.arrow.width=.5, edge.color="light gray", 
     vertex.size=20, main="Walk Trap: 4 Steps", margin=0)

layout=layout.fruchterman.reingold(edge_graph2) 
plot(edge_graph2, layout=layout, #note the use of layout
     vertex.color=cluster2_c, edge.arrow.size=.5,
     edge.arrow.width=.5, edge.color="light gray", 
     vertex.size=20, main="Walk Trap: 4 Steps", margin=0)

layout=layout.fruchterman.reingold(edge_graph_sem) 
plot(edge_graph_sem, layout=layout, #note the use of layout
     vertex.color=cluster2_sem_c, edge.arrow.size=.5,
     edge.arrow.width=.5, edge.color="light gray", 
     vertex.size=20, main="Walk Trap: 4 Steps", margin=0)
plot(edge_graph_sem, layout=layout, #note the use of layout
     vertex.color=cluster2_sem_c2, edge.arrow.size=.5,
     edge.arrow.width=.5, edge.color="light gray", 
     vertex.size=20, main="Walk Trap: 3 Steps", margin=0)


edge_graph_eb=cluster_edge_betweenness(graph=edge_graph)
edge_graph_eb_sem=cluster_edge_betweenness(graph=edge_graph_sem)

mems_eb=membership(edge_graph_eb)
mems_eb_sem=membership(edge_graph_eb_sem)

# LinkRank 
lr.modularity <- function(g,
                          partition, 
                          damping = .85, 
                          pr.algo = 'prpack',
                          weights = NULL) {
  
  ## g           = graph (igraph object)
  ## partition   = graph partition (numeric vector of memberships or "communities" object)
  ## damping     = damping factor (1 - teleportation prob.)
  ## pr.algo     = algorithm to calculate Perron vector,
  ##               possible options are "prpack", "arpack", and "power"
  ## weights     = If this is NULL and the graph has a weight edge attribute
  ##               then that is used. If weights is a numerical vector then
  ##               it used, even if the graph has a weights edge attribute.
  ##               If this is NA, then no edge weights are used (even if the
  ##               graph has a weight edge attribute)
  
  # check args
  if (!is.igraph(g)) 
    stop('graph is not an i.graph object')
  
  if (damping > 1 | damping < 0) 
    stop('damping factor has to be between zero and one!')
  
  # get algorithm name to calculate Perron vector
  pr.algo <- match.arg(pr.algo, c('prpack','arpack','power'))
  
  # no of nodes
  n <- vcount(g)
  # node sequence
  v.seq <- seq_len(n)
  
  # get membership vector
  if (class(partition) == 'communities') {
    
    pp <- membership(partition)
    
  } else {
    
    if (!is.numeric(partition))
      stop("'partition' has to be a 'communities' object or a numeric vector!")
    pp <- partition
    
  }
  
  # check dimensions
  if (length(pp) != n) 
    stop('Length of membership vector differs from number of nodes!')
  
  # get adjacency matrix & out-degree
  if (is.vector(weights) & length(weights) > 1) {
    
    # check args
    if (ecount(g) != length(weights))
      stop("'weights' differes in length from ecount!")
    if (!is.numeric(weights))
      stop("'weights' must be 'NA','NULL', or a numeric vector!")
    
    edge_attr(g, 'tmp') <- weights
    A <- get.adjacency(g, type = 'both', attr = 'tmp')
    
    out.deg <- strength(g, mode = 'out', weights = weights)
    
  } else if (is.null(weights)) {
    
    if ('weight' %in% edge_attr_names(g)) {
      
      A <- get.adjacency(g, type='both', attr='weight')
      out.deg <- strength(g, mode = 'out')
      
    }  else {
      
      A <- get.adjacency(g, type='both')
      out.deg <- degree(g, mode = 'out')
      
    }
    
  } else if (is.na(weights)) {
    
    A <- get.adjacency(g, type='both')
    out.deg <- degree(g, mode = 'out')
    
  } else {
    
    stop("'weights' option has to be 'NA','NULL', or a numeric vector!")
    
  }
  
  # dead-end nodes
  dangling <- out.deg == 0
  
  # row-normalize A (recycle vector)
  G.temp <- A / out.deg
  # equivalent to sweep(A, 1, out.deg, FUN='/')
  
  # set rows for dead-end nodes to zero
  if (sum(dangling) > 0) {
    G.temp[dangling,] <- 0
  }
  
  # add teleportation probabilities
  Tmat <- Matrix::Matrix(1/n * (damping * dangling + 1 - damping), 
                         nrow = n, ncol = n)
  G <- damping * G.temp + Tmat
  
  # get Perron vector (PageRank)
  p.vec <- page_rank(g, damping = damping, algo = pr.algo, weights = weights)$vector
  
  # LinkRank matrix
  Q <- G * p.vec -  tcrossprod(p.vec)
  # equivalent to sweep(G, 1, p.vec, '*') -  tcrossprod(p.vec)
  
  # get LinkRank Modularity by summing over within-community weights
  return(sum(Q[outer(pp, pp, '==')]))
  
}

edge_graph_un=graph_from_data_frame(d=sem_edge_di1, directed=F,
                                 vertices=sem_attr1) 


infomap_1 <- cluster_infomap(edge_graph_un, e.weights = NULL, v.weights = NULL,
                nb.trials = 10, modularity = TRUE)
membership(infomap_1)
communities(infomap_1)

layout=layout.fruchterman.reingold(edge_graph_un) 
plot(edge_graph_un, layout=layout, #note the use of layout
     vertex.color=infomap_1, edge.arrow.size=.5,
     edge.arrow.width=.5, edge.color="light gray", 
     vertex.size=20, main="Infomap", margin=0)

