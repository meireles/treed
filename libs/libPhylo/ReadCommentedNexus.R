# This function parses a Nexus file and keeps the comments
# It mostly uses RBrownie's function read.simmap.new() after modifying the comment structure

require(phyext,quietly=T)
require(phylobase,quietly=T)
require(ade4, quietly=T)

readCommentedNexusTree <- function(file){
  # Read in the text file
  x <- readLines(file, warn = F)
  
  # Clean-up tabs
  x <- gsub("\t","",x)
  
  # Figure out where the tree block is and cut the rest off
  treeblockStart <- grep("begin trees;",x,ignore.case=T) + 1 # the + 1 skips "begin trees;" itself
  treeblockEnd <- grep("end;",x, ignore.case=T)
  treeblockEnd <- treeblockEnd[ which(treeblockEnd>treeblockStart)[1]  ] - 1  # skips "end;" itself
  x <- x[ seq(treeblockStart, treeblockEnd) ]
  
  # Get Tip labels in a translation table if there is one
  getTranslTable <- function(x){
    start <- which(regexpr("Translate",x) >=1 )[ 1 ]
    end <- which(regexpr(";",x)>=1)[1]
    translation <- x[c((start+1):(end-1))]
    translation <- gsub(",","",translation)
    translation
  }
  
  taxonNames <- getTranslTable(x)
  
  # Find the tree statments, store and count them
  treeStatement <- x[grep("tree",x,ignore.case=T)]
  nTrees <- length(treeStatement)
  
  # Get the commented parenthetical statement (no tree name, [&R], etc.)
  # goes from the 1st parenthesis to the semicolon 
  treeStatement <- substr(treeStatement,min(unlist(gregexpr("\\(",treeStatement))),max(unlist(gregexpr(";",treeStatement))))
  
  # read.simmap.new likes the tree to be (A:[comment]brlens) while beast does (A[comment]:brlens)
  # just change the colon position in the tree statement	
  treeStatement <- gsub("]:","]", treeStatement)
  treeStatement <- gsub("\\[",":[", treeStatement)
  
  # Get the trees with read.simmap.new and store them in a list
  trees <- as.list(c(1:nTrees))
  for(i in 1:nTrees){
    trees[[i]] <- phyext::read.simmap.new(text=treeStatement[i])
    phylobase::tipLabels(trees[[i]]) <- taxonNames[as.numeric(phylobase::tipLabels(trees[[i]]))]
  }
  
  if(nTrees==1){
    return(trees[[1]]) # So it retunns the tree unlisted and with the right class
  } else{
    warn("This object is a simple list of phylo4 objetcs")
    return(trees) # This will return a list of trees, which has no class as a whole. Lookup a "multiphylo4" class in Phylobase.
  }
}
