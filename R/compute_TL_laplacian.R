#' Compute topological ordering
#' 
#' Method to compute trophic levels from an igraph object G with the method described in MacKay et al 2020.
#' @param G The metaweb, it has to be an igraph object
#' @return The trophic level of each node of G. By ordering species according to their trophic level, we obtain one topological ordering of the graph.
#' @references MacKay, R. S., Johnson, S., & Sansom, B. (2020). How directed is a directed network?. Royal Society open science, 7(9), 201138.
#' @importFrom igraph get.adjacency degree V vcount
#' @importFrom stats setNames
#' @importFrom Matrix solve
#' @author Giovanni Poggiato
#' @examples 
#' data(G)
#' compute_TL_laplacian(G)
#' @export
compute_TL_laplacian <- function(G){
  # recursive function on connected components of G
  if (igraph::vcount(G) == 1) return(setNames(0, igraph::V(G)$name))
  A = as.matrix(igraph::get.adjacency(G))
  names_loc = rownames(A)
  
  
  if (isNamespaceLoaded("sna")) { try(unloadNamespace("sna"), silent = TRUE) }
  if (!isNamespaceLoaded("igraph")) { 
    if(!requireNamespace('igraph', quietly = TRUE)) stop("Package 'igraph' not found")
  }
  
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



