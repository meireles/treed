### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### %%%
### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xyPhyloCoord <- function(phy){
  if(class(phy)!="phylo"){
    phy <- as(phy,"phylo")
  }
  # Grab some basic info about the phylogeny
  ntip <- Ntip(phy)
  nnode <- Nnode(phy,internal.only=TRUE)
  totnode <- Nnode(phy, internal.only=FALSE)
  root <- ntip + 1
  # Set the XX coords to be the distance between each node and the root
  xx <- dist.nodes(phy)[,root]
  # Calculate the YY coordinates
  yy <- rep(NA,totnode)
  cur <- 0
  for(i in rev(root:totnode)){
    desc <- phy$edge[which(phy$edge[,1]==i),2]
    tips <- desc[which(desc<=ntip)]
    if(length(tips) > 0){
      cur <- seq(length(tips))+max(cur)
      yy[tips] <- cur
    }		
    yy[i] <- mean(yy[desc])
  }
  xy <- cbind(xx=xx,yy=yy/ntip)
  return(xy)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### %%%
### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotNodes3d <- function(xyz,radius=.01,col="black",...){
  rgl.spheres(x=xyz, radius=radius, col=col,...)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### %%%
### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotEdges3d <- function(phy,xyz,style=c("square","diagonal"),edgeType=c("line","cylinder"),
                        cylSides=16,cylRadius=.02,lineLwd=1,col="black",...){
  
  if(class(phy)!="phylo"){
    phy <- as(phy,"phylo")
  }
  
  if(length(col)==1){
    col <- rep(col,nrow(xyz))
  }
  
  # Get the edge matrix. Pretty unecessary with the newer code. Take care of it!
  edges <- phy$edge
  # Generate a list of edges with makeEdges3d..
  edgeList <- makeEdges3d(phy=phy,xyz=xyz,style=style)
  
  # Plotting these edges in a loop is bad. Normally in opengl I would generate a list 
  # of edges and then call the plotting functions only once, but I cannot get that
  # to work well with rgl...
  if(style == "diagonal"){
    if(edgeType == "line"){
      for(i in 1:nrow(edges)){
        rgl.lines(edgeList[[i]],lwd=lineLwd,col=col[i],...)
      }	
    }
    if(edgeType == "cylinder"){
      for(i in 1:nrow(edges)){
        shade3d(cylinder3d(edgeList[[i]], radius=cylRadius,sides=cylSides,
                           e2=rbind(c(1,1,1),c(1,1,1))), col=col[i])
      }	
    }
  }
  if(style == "square"){
    if(edgeType == "line"){
      for(i in 1:nrow(edges)){
        rgl.lines(edgeList[[i]]$vBranch,lwd= lineLwd,col=col[i],...)
        rgl.lines(edgeList[[i]]$hBranch,lwd= lineLwd,col=col[i],...)
      }
    }
    if(edgeType == "cylinder"){
      for(i in 1:nrow(edges)){
        v <- cylinder3d(edgeList[[i]]$vBranch, radius=cylRadius,sides= cylSides,
                        e2=rbind(c(1,1,1),c(1,1,1)))
        h <- cylinder3d(edgeList[[i]]$hBranch, radius=cylRadius,sides= cylSides,
                        e2=rbind(c(1,1,1),c(1,1,1)))
        # Add normals to smooth the objects...
        v <- addNormals(v)
        h <- addNormals(h)
        # And plot them
        shade3d(v, col=col[i])
        shade3d(h, col=col[i])
        # The following code is very bad design. I now that vBranch[2,] = hBranch[1,] = the "elbow" on a square edge
        # I think that I should define the "elbow" node somewere instead of just calling it "indirectlly. Anyways, here
        # i'll plot spheres on the elbows with the same color and diameter of the cylinders...
        material3d(smooth=TRUE,specular="black")
        rgl.spheres(x= edgeList[[i]]$vBranch[2,], radius=cylRadius, col=col[i],...)
      }
    }
  }
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### %%%
### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
makeEdges3d <- function(phy,xyz,style=c("square","diagonal"),...){
  
  if(class(phy)!="phylo"){
    phy <- as(phy,"phylo")
  }
  
  edges <- phy$edge
  
  ### TEST
  edges <- phy$edge[order(phy$edge[,2]),]
  
  
  a <- 1	#ancestor column on edges matrix
  d <- 2	#descendent column on edges matrix
  tmp <- list()
  edgeList<- list()
  # This is a bad loop...			
  if(style == "diagonal"){
    for(i in 1:nrow(edges)){
      tmp <- rbind(xyz[edges[i,1],], xyz[edges[i,2],])
      edgeList[[i]] <- tmp
    }
  }
  if(style == "square"){
    for(i in 1:nrow(edges)){
      # Make to branches, a "horizontal" and a "vertical".
      # Calculate the "elbow nodes"
      vBranch <- cbind(
        c(xyz[edges[i,a],"xx"],xyz[edges[i,a],"xx"]), 
        c(xyz[edges[i,a],"yy"],xyz[edges[i,d],"yy"]), 
        c(xyz[edges[i,a],"zz"],xyz[edges[i,a],"zz"]))
      hBranch <- cbind(
        c(xyz[edges[i,a],"xx"],xyz[edges[i,d],"xx"]),
        c(xyz[edges[i,d],"yy"],xyz[edges[i,d],"yy"]),
        c(xyz[edges[i,a],"zz"],xyz[edges[i,d],"zz"]))
      tmp <- list(vBranch= vBranch,hBranch= hBranch)
      edgeList[[i]] <- tmp
    }
  }
  return(edgeList)
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### %%%
### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

plotLabels3d <- function(phy,xyz,tip.only=TRUE,color="black",cex=1,adj=0,...){
  
  if(class(phy)!="phylo"){
    phy <- as(phy,"phylo")
  }
  
  if(tip.only==FALSE){
    rgl::text3d(xyz,text=c(phy$tip.label, phy$node.label),...)
  }
  else{
    rgl::text3d(xyz[1:Ntip(phy),],text=phy$tip.label,color=color,cex=cex,adj=adj,...)
  }
}

### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
### %%%
### %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

printPhylo3d <- function(
  phy, xyz, fileName="name.wrl",
  setOutDir = "~/Desktop/",
  scale,
  style = c("square","diagonal"),
  nodeRadius=0.1, nodeCol="red",
  edgeLwd, edgeCol="black",...
  ){
  
  if(class(phy)!="phylo"){
    phy <- as(phy,"phylo")
  }
  
  # Set working directory and open the vrml device
  setwd(setOutDir)
  vrml.open(file = fileName, scale= scale, navigation = "PAN")
  # Get the edge matrix
  edges <- phy$edge
  # Generate a list of edges with makeEdges3d..
  edgeList <- makeEdges3d(phy=phy,xyz=xyz,style=style)
  # Print Nodes to the vrml device
  vrmlgen::points3d(x=xyz, scale=nodeRadius, col=nodeCol)
  # Print Edges to the vrml device
  if(style == "diagonal"){
    for(i in 1:nrow(edges)){
      vrmlgen::lines3d(edgeList[[i]],lwd=edgeLwd,col=edgeCol,...)
    }
  }
  if(style == "square"){
    elbows <- matrix(ncol=3,nrow=length(edgeList))
    for(i in 1:nrow(edges)){
      vrmlgen::lines3d(edgeList[[i]]$vBranch,lwd=edgeLwd,col=edgeCol,...)
      vrmlgen::lines3d(edgeList[[i]]$hBranch,lwd=edgeLwd,col=edgeCol,...)
      elbows[i,] = edgeList[[i]]$vBranch[2,]
    }
    vrmlgen::points3d(x=elbows, scale= edgeLwd,col=edgeCol,...)
  }
  vrml.close()
}