######################################################

require(Rcpp)

compute_cluster_hierarchy=function(infile, method = "STANDARD")
{
  ftable<-read.csv(infile)

  group1<- colnames(ftable)[2:ncol(ftable)]
  group2<- as.character(ftable[,1])

  count_edges <- 0
  edge_weight <- 0

  weighting_scheme <-"UNWEIGHTED"

  for( i in 1:nrow(ftable))
    for( j in 2:ncol(ftable))
      if( as.numeric(ftable[i,j]) != 0)
      {
        count_edges <- count_edges+1
      
        if(count_edges == 1)
        {
          edge_weight <- as.numeric(ftable[i,j])
        } else {
          
         if(edge_weight != as.numeric(ftable[i,j]))
           weighting_scheme <- "WEIGHTED" 
        }

      } # if( as.numeric(ftable[i,j]) != 0)

  flags <- list(weighting_scheme,method)

  edge_nodes1 <-vector(mode="integer",count_edges)
  edge_nodes2 <-vector(mode="integer",count_edges)
  edge_weights <-vector(mode="numeric",count_edges)

  count_edges <- 1

  for( i in 1:nrow(ftable))
    for( j in 2:ncol(ftable))
      if( as.numeric(ftable[i,j]) != 0)
      {
        edge_nodes1[count_edges] <- j-1
        edge_nodes2[count_edges] <- i
        edge_weights[count_edges] <- as.numeric(ftable[i,j])  
        count_edges<-count_edges+1
      }

  count_edges<-count_edges-1

  n1 <- length(group1)
  n2 <- length(group2)
  cluster_graph_x <-vector(mode="integer",count_edges)
  cluster_graph_y <-vector(mode="numeric",count_edges)
  ordered_edges <-vector(mode="integer",count_edges)
  parent_nodes <-vector(mode="integer",(2*count_edges)-1)
  depths <-vector(mode="integer",(2*count_edges)-1)
  edge_ranges_begin <-vector(mode="integer",(2*count_edges)-1)
  edge_ranges_end <-vector(mode="integer",(2*count_edges)-1)
  number_of_cluster_nodes <- -1

  for( i in 1:count_edges)
  {
    cluster_graph_x[i] <- -1
    cluster_graph_y[i] <- -1
  }

  sourceCpp("construct_community_hierarchy.cpp")

  res <-construct_community_hierarchy(as.character(flags), as.character(group1), as.character(group2), as.integer(edge_nodes1), as.integer(edge_nodes2), as.numeric(edge_weights) )

  ordered.edges <- res[[1]]

  new_edge_nodes1 <-vector(mode="integer",count_edges)
  new_edge_nodes2 <-vector(mode="integer",count_edges)
  new_edge_weights <-vector(mode="numeric",count_edges)

  for(i in 1:count_edges)
  {
    new_edge_nodes1[i] = edge_nodes1[ordered.edges[i]] 
    new_edge_nodes2[i] = edge_nodes2[ordered.edges[i]]
    new_edge_weights[i] = edge_weights[ordered.edges[i]]
  } 

  cluster_tree <- list( node.count.group1 = n1, node.group1 = group1, node.count.group2 = n2, node.group2 = group2, edge.count = length(new_edge_nodes1), edge.nodes1 = new_edge_nodes1, edge.nodes2 = new_edge_nodes2, parent.nodes = res[[2]], node.depths = res[[3]], edge.range.begin = res[[4]], edge.range.end = res[[5]], number.of.cluster.nodes = length(res[[2]]), cluster.graph.x = res[[6]], cluster.graph.y = res[[7]])

  ###############################################################################
  # Output the graph with the performance values for each partition of clusters #
  ###############################################################################

  range_index <- 1

  cgg_x <- cluster_tree[['cluster.graph.x']]
  cgg_y <- cluster_tree[['cluster.graph.y']]

  if(length(cgg_x) == 0)
  {
    print("The input dataset contained no network nodes")
  } else {

    if( method == "STANDARD" )
    {
      explanation <- "Similarity function used: Tanimoto Coefficient"
    } else {
    if( method == "MIN_WEIGHT_SUM" )
      explanation <- "Similarity function used: Minimum Weight Sum"
    }

    plot(cgg_x,cgg_y,"b",xlab="Number of Clusters",ylab="Partition Density", main = explanation)
  } 

  return(cluster_tree)

} #compute_cluster_hierarchy=function(...)

extract_clusters=function(cluster_tree,depth)
{
  parents <- cluster_tree[['parent.nodes']]
  depths <- cluster_tree[['node.depths']]
  edge_range_begin <- cluster_tree[['edge.range.begin']]
  edge_range_end <- cluster_tree[['edge.range.end']]
  edge_nodes1 <- cluster_tree[['edge.nodes1']]
  edge_nodes2 <- cluster_tree[['edge.nodes2']]
  node_group1 <- cluster_tree[['node.group1']]
  node_group2 <- cluster_tree[['node.group2']]

  number_of_clusters <-0

  # Due to impractilities on building lists in R incrementally,
  # we first count all the elements that are to be stored in the 
  # triple list of the output   

  for( i in 1:cluster_tree[['number.of.cluster.nodes']])
    if( depths[i] >= depth && ( parents[i]==0 || depths[parents[i]] < depth) )
      number_of_clusters <- number_of_clusters + 1

  temp_count <- 0
  groups1<-vector(mode="integer",number_of_clusters)
  groups2<-vector(mode="integer",number_of_clusters)

  for( i in 1:cluster_tree[['number.of.cluster.nodes']])
    if( depths[i] >= depth && ( parents[i]==0 || depths[parents[i]] < depth) )
    {
      temp_count <- temp_count+1  
      groups1[temp_count] = length(unique(edge_nodes1[edge_range_begin[i]:edge_range_end[i]]))
      groups2[temp_count] = length(unique(edge_nodes2[edge_range_begin[i]:edge_range_end[i]]))
    }

  clusters <-list(list(list()))

  temp_count = 0

  for( i in 1:cluster_tree[['number.of.cluster.nodes']])
  {
    if( depths[i] >= depth && ( parents[i]==0 || depths[parents[i]] < depth) )
    {
      temp_count <- temp_count + 1
      subgroup1 <- edge_nodes1[edge_range_begin[i]:edge_range_end[i]]
      subgroup2 <- edge_nodes2[edge_range_begin[i]:edge_range_end[i]]

      sg1 <-unique(subgroup1)
      sg2 <-unique(subgroup2)

      subg1 <-order(sg1) 
      subg2 <-order(sg2)

      species1 <-vector("character", groups1[temp_count])
      species2 <-vector("character", groups2[temp_count])

      species1 <-node_group1[sg1[subg1]]
      species2 <-node_group2[sg2[subg2]]

      groups <-list(species1,species2)
      names(groups)[1] <- "Group 1"
      names(groups)[2] <- "Group 2"

      clusters[[temp_count]] <- groups
      names(clusters)[temp_count] = paste("Cluster",temp_count,sep = " ")
    }
 
  } #for( i in 1:cluster_tree[['number.of.cluster.nodes']])

  return(clusters)

} #extract_clusters=function(...)


# Methods: "STANDARD" or "MIN_WEIGHT_SUM" 

##################
# Example Usage ##
##################

#infile<-"~/Desktop/Species_network_R_interface/Data_sets/SorenBek.csv"
#cutdepth <- 10
#ct<-compute_cluster_hierarchy(infile,"STANDARD")
#original<-extract_clusters(ct,cutdepth)
#print(original)

