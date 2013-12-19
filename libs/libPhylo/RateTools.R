
getAverageMolEvoRate <- function(phy, rate){

  #' Computes the average rate of evolution weighted by branch lengths
  #'
  #' This function computes the an average rate of evolution weigted by branch 
  #' lengths.
  #' 
  #' @param phy a single tree in ape or phylobase formats (that is, of classes \
  #' phylo or phylo4/d object 
  #' @param rate a vector specifying the rates of rates of evolution for every \
  #' edge in the phylogeny. This vector must be in the same order as the phylogeny's \
  #' edge matrix, i.e. phy$edge (ape) or phy@edge (phylobase)
  #' @return computed weighted evolution rate (float)
  #'
  #' @examples
  #' #TODO
  
  # First check if all the branches have rates
  if(length(rate[is.na(rate)]) >= 2){
    stop("Not all branches have rates! Only the root can be an NA \n Check your 
         rate vector or your treeAnnotator settings")
  }
  
  if(class(phy) == "phylo") {	
    stop("not implemented yet")
    # suppressPackageStartupMessages(require(ape,quietly=T))
    # return(
    # sum(rate * phy$edge.length) / sum(phy$edge.length)
    # )
  }
 
  if(class(phy) == "phylo4" | class(phy) == "phylo4d"){
    suppressPackageStartupMessages( require(phylobase, quietly=T) )
    brlens <- orderBrLensByNodeName(phy)
    return(
      sum(na.omit(rate * brlens)) / sum(na.omit(brlens))
      )
  }
  else{
    stop("The phy object provided is not of any known tree class")
  }	
}

getDeltaMolEvoRate <- function(phy, rate){
  # Calulates the delta molecular evolution rate as the difference between each 
  # branch's rate and the treewide avg rate.
  # Takes in:
  # phylo or phylo4/d object
  # and a vector of rates of evolution
  # RATES MUST be in the SAME ORDER as they appear in the EDGE MATRIX
  # i.e. phy#edge (ape) or phy@edge (phylobase)
  # see "getDataFromPhylo4d" to see how to reorder data from a phylo4d
  # Returns:
  # vector of delta rates of evolution; w/ length = number of edges
  return(
    rate - getAverageMolEvoRate(phy, rate)
    )
}

getAccumSubsAtNode <- function(phy,rate){
  # Calulates the acumulated number of substitutions at a given node
  # Takes in:
  # phylo or phylo4/d object
  # and a vector of rates of evolution, which can be:
  # untransformed: straight from BEAST, e.g. "phy@data$rate"
  # transformed to be a delta rate, e.g. "getDeltaMolEvoRate()"
  # RATES MUST be in the SAME ORDER as they appear in the EDGE MATRIX
  # i.e. phy#edge (ape) or phy@edge (phylobase)
  # see "getDataFromPhylo4d" to see how to reorder data from a phylo4d
  # Returns: 
  # vector of accumulated substitutions at a given node
  
  if(class(phy) == "phylo4" | class(phy) == "phylo4d"){
    suppressPackageStartupMessages(require(phylobase,quietly=T))
    
    brlens <- orderBrLensByNodeName(phy)
    
    subsPerEdge <- brlens * rate
    subsPerEdge[is.na(subsPerEdge)] <- 0
    
    anc <- phylobase::ancestors(phy,phy@edge[,"descendant"],type="ALL")
    
    accum <- rep(NA,phylobase::nEdges(phy))
    # For the next for loop, the accum vector will have the same order as anc and edge matrix
    # I'll reorder it as increasing node number before output
    names(accum) <- as.numeric(names(anc))
    for(i in 1:phylobase::nEdges(phy)){
      accum[i] <- sum(subsPerEdge[anc[[i]]])
    }
    accum <- accum[order(as.numeric(names(accum)))] # Reordering (see comment above)
    return(accum)
  }
}


orderBrLensByNodeName <- function(phy){
  br <- phy@edge.length
  names(br) <- as.numeric(phy@edge[,"descendant"])
  br <- br[order(as.numeric(names(br)))]
  return(br)
}
