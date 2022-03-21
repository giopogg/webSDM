##########################################################################################################
#### Play with the metaweb to find a dag
#### G. Poggiato 18/02/2021
library(tidyr)
library(dplyr)
library(igraph)
library(NetIndices)
library(plyr)

rm(list=ls())
dir="/Users/poggiatg/Documents/Phd/Futureweb/Code/"
setwd(dir)

metaweb=read.csv2(paste0(dir,"DATA/metaweb/TetraEU_pairwise_interactions.csv"))
SppNames=read.csv2(paste0(dir,"DATA/metaweb/Id_SppName.csv"))

metaweb = dplyr::select(metaweb,sourceTaxonId,sourceLifestageName,targetTaxonId,targetLifestageName)

#remove cannibalism
metaweb = dplyr::filter(metaweb, !sourceTaxonId==targetTaxonId)


#remove non adults-adults interactions
metawebAdult = dplyr::filter(metaweb,targetLifestageName=="adults")

summarise(metawebAdult,length(targetLifestageName))

#build graph for adult-adult only

MetaGraph = graph_from_edgelist(as.matrix(dplyr::select(metawebAdult,sourceTaxonId,targetTaxonId)), 
                                directed = TRUE)

MetaAdj = as_adjacency_matrix(MetaGraph)

save(MetaAdj,file="DATA/metaweb/MetaAdjAdult.RData")

#TO SAVE AND OPEN IN GEPHI

to_be_saved = as.data.frame(as_edgelist(MetaGraph, names = T))
colnames(to_be_saved)=c("Source","Target")
rownames(to_be_saved)=NULL
write.csv(to_be_saved,file="DATA/metaweb/GephiMetaAdult.csv", row.names = F)

SppNamesAdult = dplyr::filter(SppNames, Id %in% metawebAdult$sourceTaxonId | Id %in% metawebAdult$targetTaxonId )
write.csv(SppNamesAdult,file="DATA/metaweb/SppNamesAdult.csv",row.names = F)

# The whole metaweb is not a DAG
is.dag(MetaGraph)

##########################################################################################################
#### Select Mammals only

metawebMam = filter(metawebAdult, grepl("M", sourceTaxonId), grepl("M", targetTaxonId))

MetaGraphMam = graph_from_edgelist(as.matrix(dplyr::select(metawebMam,sourceTaxonId,targetTaxonId)), 
                                directed = TRUE)
is.dag(MetaGraphMam)

plot(MetaGraphMam)

MetaAdjMam = as_adjacency_matrix(MetaGraphMam)

save(MetaAdjMam,file="DATA/metaweb/MetaAdjMamAdult.RData")

#TO SAVE AND OPEN IN GEPHI
to_be_saved = as.data.frame(as_edgelist(MetaGraphMam, names = T))
colnames(to_be_saved)=c("Source","Target")
rownames(to_be_saved)=NULL
write.csv(to_be_saved,file="DATA/metaweb/GephiMetaMamAdult.csv", row.names = F)

SppNamesMam = filter(SppNames, Id %in% metawebMam$sourceTaxonId | Id %in% metawebMam$targetTaxonId )
write.csv(SppNamesMam,file="DATA/metaweb/SppNamesMam.csv",row.names = F)

TLMam=TrophInd(as.matrix(MetaAdjMam))

hist(TLMam$TL)

#Max TL is Tadarida_teniotis
SppNames[SppNames$Id==rownames(filter(TLMam, TL==max(TLMam$TL))),]
#Max OI is Lepus_capensis
SppNames[SppNames$Id==rownames(filter(TLMam, OI==max(TLMam$OI))),]

# Some topological ordering

SortedMam = topo_sort(MetaGraphMam,mode="in")
SortedMam.df = data.frame(Id=names(SortedMam))
join(SortedMam.df,SppNames,by="Id")

SppNames[SppNames$Id %in% names(SortedMam),]
SppNames[SppNamesAdult$Id %in% preys$targetTaxonId,]

#How many Spp are basal?
basal = filter(SppNamesMam, !(Id %in% metawebMam$sourceTaxonId))
basal = SppNamesMam[which(SppNames$Id %in% metawebMam$sourceTaxonId),]
dim(SppNamesMam)
dim(basal) #only 38 are not basal ahah

filter(metaweb, sourceTaxonId =="M101")
##########################################################################################################
#### Louise's predators subset

# top_pred=c("B88","B103","B111","B123","R106","R79","R208","M80","M31","M165","M111")
# 
# # Finds all vertices reachable from a given vertex (when mode="out'), or the opposite: all vertices from
# #  which a givenvertex is reachable via a directed path (when mode = 'in'). When mode="out" it assumes
# #  a directed graph (all= indirected graph)
# 
# all_preys = subcomponent(MetaGraph,top_pred,mode="out")
# 
# induced_from_preys = induced_subgraph(MetaGraph,all_preys)
# 
# is.dag(induced_from_preys)
# 
# preysNb = metawebAdult %>% filter(sourceTaxonId %in% top_pred) %>% dplyr::select(targetTaxonId) %>% unique
# length(preysNb$targetTaxonId)
# #still 1018 spp...
# 
# #reptiles (and lupus?) way too generalist
# degree(MetaGraph)[top_pred]
# 
# #even with just mammals it's not enough
# top_pred=c("M80","M165","M111")
# 
# # Finds all vertices reachable from a given vertex, or the opposite: all vertices from which a given 
# # vertex is reachable via a directed path. When mode="out" it assumes a directed graph 
# # (all= indirected, in= directed but with inverted arrows.)
# 
# all_preys = subcomponent(MetaGraph,top_pred,mode="out")
# induced_from_preys = induced_subgraph(MetaGraph,all_preys)
# 
# length(V(induced_from_preys))
# 
# is.dag(induced_from_preys)


#######################################################################################################
##### Find an easy toy dataset


#top_pred=c("B88","B103","B111","B123","R106","R79","R208","M80","M31","M165","M111")


#start from a few interesting species and select some of their interesting preys
# preys = metawebAdult %>% filter(sourceTaxonId=="B103") %>% dplyr::select(targetTaxonId)
# 
# SppNamesAdult[SppNamesAdult$Id %in% preys$targetTaxonId,]

#B287             Dendrocopos_major
# M128             Microtus_agrestis

#It's a mess, probably better to jjust take mammals


#############################################################################################################
############ Why it's non dag?



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


MetaGraph_TL = compute_TL_laplacian(MetaGraph)
hist(MetaGraph_TL)
sortSp = names(sort(MetaGraph_TL,decreasing=T))

#Build a complete df of problematic links
notDAG_vec=vector()
notDAG=setNames(data.frame(matrix(ncol = 6, nrow = 0)), c("Pred", "Prey", "Direct_loop","Three_sp_loop","Less_than_8_sp_loop","High_dim_loop"))

for(i in sortSp){
  #what i feeds on that have a higher TL than i
  notDAG_i = names(which(MetaGraph_TL[neighbors(MetaGraph,i,mode=c("out"))] < MetaGraph_TL[i]))
  i_Sp = SppNamesAdult$Label[which(SppNamesAdult$Id==i,)]
 
   if(length(notDAG_i>0)){
    for(j in notDAG_i){
    
      ij_Sp = SppNamesAdult$Label[which(SppNamesAdult$Id == j,)]
      
        if(MetaAdj[j,i]==1){ #if j eats i, than i<->j

              notDAG=rbind(notDAG,data.frame(Pred=i_Sp,Prey=ij_Sp, Direct_loop=T,Three_sp_loop=F,Four_sp_loop=F,High_dim_loop=F))
      
        }else{#check if 3spp loop
                loop3dSpp=intersect(names(neighbors(MetaGraph,j,mode=c("out"))),names(neighbors(MetaGraph,i,mode=c("in"))))
                if(length(loop3dSpp)>0){
                notDAG=rbind(notDAG,data.frame(Pred=i_Sp,Prey=ij_Sp, Direct_loop=F,Three_sp_loop=paste(SppNamesAdult$Label[which(SppNamesAdult$Id %in% loop3dSpp,)],collapse=","),Four_sp_loop=F,High_dim_loop=F))
                }else{
                  allPaths=all_simple_paths(MetaGraph,from=j,to=i,mode="out",cutoff=3)
                  
                  if(length(allPaths)>0){   
                    notDAG=rbind(notDAG,data.frame(Pred=i_Sp,Prey=ij_Sp, Direct_loop=F,Three_sp_loop=F,Four_sp_loop=T,High_dim_loop=F))
                  }else{
                    notDAG=rbind(notDAG,data.frame(Pred=i_Sp,Prey=ij_Sp, Direct_loop=F,Three_sp_loop=F,Four_sp_loop=F,High_dim_loop=T))
                    
                  }
                }
        }
      
      notDAG_vec=c(notDAG_vec,paste0(ij_Sp,",",i_Sp))
      
    }
    
  }
    
}

write.csv(notDAG,file=paste0(dir,"DATA/metaweb/ProbLinkDAG.csv"))


#Code to check species
i="M114"
i="B445"
i="B257"
i="B178"
#j = names(which(MetaGraph_TL[neighbors(MetaGraph,i,mode=c("out"))] < MetaGraph_TL[i]))
i="R213"
SppNamesAdult$Label[which(SppNamesAdult$Id %in% names(neighbors(MetaGraph,i,mode=c("in"))))]
SppNamesAdult$Label[which(SppNamesAdult$Id %in% names(neighbors(MetaGraph,i,mode=c("out"))))]

# i="B257" Piu animali mangiano animali che mangiano il gufo (volpe e mangusta).
# for(j in names(neighbors(MetaGraph,i,mode=c("in")))) {
#   print(paste0("xx",j,"xx"))
# print(SppNamesAdult$Label[which(SppNamesAdult$Id %in% names(neighbors(MetaGraph,j,mode=c("in"))))])
# }

# nessuno mangia qualcuno che mangia la mangusta
# i="M88"
# for(i in names(neighbors(MetaGraph,i,mode=c("in")))) print(paste0("xx",i,"xx")); print(SppNamesAdult$Label[which(SppNamesAdult$Id %in% names(neighbors(MetaGraph,i,mode=c("in"))))])

#j = names(which(MetaGraph_TL[neighbors(MetaGraph,i,mode=c("out"))] < MetaGraph_TL[i]))




#What if we remove all self loops and 3d loops first
MetaDAG_adj=MetaAdj

for(i in sortSp){
  #what i feeds on that have a higher TL than i
  notDAG_i = names(which(MetaGraph_TL[neighbors(MetaGraph,i,mode=c("out"))] < MetaGraph_TL[i]))

  if(length(notDAG_i>0)) {
    for(j in notDAG_i){
      #MetaDAG_adj[i,j]=0
      if(MetaAdj[j,i]==1){ #if j eats i, than i<->j
        MetaDAG_adj[j,i]=0

      }
    }
  }
}
is.dag(graph_from_adjacency_matrix(MetaDAG_adj))

Meta_no_sf=graph_from_adjacency_matrix(MetaDAG_adj)

# And then we recheck
MetaGraph_TL = compute_TL_laplacian(Meta_no_sf)
hist(MetaGraph_TL)
sortSp = names(sort(MetaGraph_TL,decreasing=T))

notDAG=vector()
for(i in sortSp){
  #what i feeds on that have a higher TL than i
  notDAG_i = names(which(MetaGraph_TL[neighbors(Meta_no_sf,i,mode=c("out"))] < MetaGraph_TL[i]))

  if(length(notDAG_i)>0) {
  notDAG_i_Sp = SppNamesAdult$Label[which(SppNamesAdult$Id %in% notDAG_i,)]
  i_Sp = SppNamesAdult$Label[which(SppNamesAdult$Id==i,)]


  notDAG=c(notDAG,paste0(notDAG_i_Sp,",",i_Sp))
  }

}
notDAG
