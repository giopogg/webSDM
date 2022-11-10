#' Builds SDM formulae
#'
#' Builds the formula of both the abiotic and biotic terms to fit a single species SDM based on the input parameters. The function is called inside the SDMfit function
#' @param form.init The abiotic part of the formula
#' @param species The preys (or predators) of the focal species
#' @param sp.formula optional parameter for composite variables. See ?trophicSDM
#' @param sp.partition optional parameter to specify groups of species for composite variables. See ?trophicSDM
#' @param useBRMS whether brms is used (TRUE if penal = "coeff.signs" and method = "stan_glm).
#' @author Giovanni Poggiato and Jérémy Andréoletti
#' @importFrom stats formula as.formula

buildFormula <- function(form.init, species, sp.formula=NULL, sp.partition=NULL, useBRMS){
  #if no composite variables are used
  if(is.null(sp.formula)){
    if(useBRMS){
      # "tempYi" are intermediary variables allowing to set lower/upper-bounds on priors
      form.all = paste(as.character(form.init), paste(paste0("temp",species), collapse="+"), sep="+")
      # intermediary species variables are defined hierarchically from observed species variables
      form.brms = lapply(species, function(sp) as.formula(paste0("temp",sp," ~ 0 + ",sp)))
    } else form.all = paste(as.character(form.init), paste(species,collapse="+"), sep="+")
  }else{
    # if no group definition
    if(grepl("richness",sp.formula)){
      if(is.null(sp.partition)){
        # we replace "richness" with the sum of all species (and all eventual other terms), both in linear and eventually other terms
        form.all=as.character(form.init)

        # replace richness by the observed species variables
        if(useBRMS){
          # hierarchical definition
          form.rich = "richness"
          form.brms = paste0("richness ~ 0 + ",
                             gsub("richness",paste0("I(",paste(paste0(species),collapse = "+"),")"),
                                  sp.formula))
        } else form.rich = gsub("richness", paste0("I(",paste(species,collapse = "+"),")"), sp.formula)
        form.all = paste(form.all,form.rich,sep="+")

        
        
        
      }else{
        # all as above but repeated for every cluster
        form.all=as.character(form.init)


        form.rich = NULL
        for(j in 1:length(sp.partition)){
          species.j = species[species %in% sp.partition[[j]]]

          if(length(species.j)>0){        # at least 2 species in the cluster
            if(is.null(form.rich)){
              # create form.rich
              if(useBRMS){
                form.rich = paste0("richness",j)
                form.brms = list(paste0("richness", j,
                                        " ~ 0 + ", gsub("richness",
                                                        paste0("(",paste(species.j,collapse="+"),")"),
                                                        paste0("I(",gsub("\\+", ")\\+I(",sp.formula),")"))))

              } else form.rich = gsub("richness",
                                      paste0("(",paste(species.j,collapse="+"),")"),
                                      paste0("I(",gsub("\\+", ")\\+I(",sp.formula),")"))
            }else{
              # add the new formula to the already existing form.rich
              if(useBRMS){
                form.rich = paste(form.rich, paste0("richness", j), sep="+")
                form.brms = c(form.brms, paste0("richness", j, " ~ 0 + ",
                                                gsub("richness", paste0("(",paste(species.j,collapse="+"),")"),
                                                     paste0("I(",gsub("\\+", ")\\+I(",sp.formula),")"))))
              } else form.rich = paste(form.rich, gsub("richness",
                                                       paste0("(",paste(species.j,collapse="+"),")"),
                                                       paste0("I(",gsub("\\+", ")\\+I(",sp.formula),")")),sep="+")
            }
          }
        }

        form.all = paste(form.all,form.rich,sep="+")

      }
    }else{

      form.all=as.character(form.init)

      #only access here if we have a list containing one formula per species. Notice that we don't check anything here, we trust the user.
      if(useBRMS){
        # hierarchical definition
        form.brms = paste0("customform ~ 0 + ", sp.formula)
        form.all = paste(form.all,"customform",sep="+")
      }else{
        form.all = as.formula(paste(form.all,sp.formula,sep="+"))
      }
    }
  }
  if (useBRMS){
    return (list(form.all=form.all, form.brms=lapply(form.brms,as.formula)))
  }
  return (list(form.all=form.all))
}
