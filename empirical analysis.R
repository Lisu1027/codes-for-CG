library(igraph)
## Function used to identify the significance of genes and obtain a graph
# estimate_sum is the estimation obtained from CG
plot_significance = function(estimate_sum,threshold1,threshold2,name,p){
  #Identify insignificant genes
  betahat = as.vector(abs(estimate_sum))
  betahat_NZ = which(betahat < threshold1)
  betahat[betahat_NZ] = 0
  estimate = matrix(betahat,p,p)
  re = estimate
  plot.adj=re
  diag(plot.adj)=0
  temp=graph.adjacency(adjmatrix=plot.adj, mode="undirected")
  temp.degree=apply(plot.adj, 2, sum)
  re1=re[temp.degree!=0,temp.degree!=0]
  
  #Identify the most significant genes
  plot.adj=re1
  diag(plot.adj)=0
  plot.adj1=plot.adj
  plot.adj1[which(plot.adj1>0)] = 1
  temp=graph.adjacency(adjmatrix=plot.adj1, mode="undirected")
  temp.degree=apply(plot.adj, 2, sum)
  V(temp)$color=(temp.degree>threshold2)+3
  length(which(V(temp)$color>3))
  name_final = name[which(temp.degree>threshold2)]
  plot(temp, vertex.size=3,
       vertex.frame.color="white",
       layout=layout.fruchterman.reingold, 
       vertex.label=NA, edge.color=grey(0.5))
  return(name_final)
}



##Function used to obtain a graph based on the results of clustering 
#Set k(the number of cluster)=3 be an example 
plot_clusters = function(p,estimate_part1,estimate_part2,estimate_part3,estimate_sum,threshold){
  betahat = as.vector(abs(estimate_sum))
  betahat_NZ = which(betahat < threshold)
  betahat[betahat_NZ] = 0
  estimate = matrix(betahat,p,p)
  re = estimate
  plot.adj=re
  diag(plot.adj)=0
  temp=graph.adjacency(adjmatrix=plot.adj, mode="undirected")
  temp.degree=apply(plot.adj, 2, sum)
  index_f = index[which(temp.degree!=0)]
  re1=re[temp.degree!=0,temp.degree!=0]
  plot.adj=re1
  diag(plot.adj)=0
  plot.adj1=plot.adj
  plot.adj1[which(plot.adj1>0)] = 1
  temp=graph.adjacency(adjmatrix=plot.adj1, mode="undirected")
  V(temp)$color=index_f
  plot(temp, vertex.size=3,
       vertex.frame.color="white",
       layout=layout.fruchterman.reingold, 
       vertex.label=NA, edge.color=grey(0.5))
}


