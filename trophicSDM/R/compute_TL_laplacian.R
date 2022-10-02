
compute_TL_laplacian <- function(G){
  # recursive function on connected components of G
  if (igraph::vcount(G) == 1) return(setNames(0, igraph::V(G)$name))
  A = as.matrix(igraph::get.adjacency(G))
  names_loc = rownames(A)
  u  = igraph::degree(G)
  v =  igraph::degree(G,mode ='in') - igraph::degree(G,mode = 'out')
  
  A = A[-1,-1]
  u = u[-1]
  v = v[-1]
  L = diag(u, ncol=length(u)) - A - t(A) # matrix is made symmetric!
  
  TL_vec = Matrix::solve(L,v)
  TL_vec = c(0,TL_vec)
  TL_vec = TL_vec - min(TL_vec)
  names(TL_vec) = igraph::V(G)$name
  return(TL_vec)
}



