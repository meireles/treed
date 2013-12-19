#' TreeD
#' Plot phylogenies in 3D using rgl and phylogenetics packages
#' @author dudu meireles
#' @email jcm54(at)duke.edu


# Load My Custom Functions
mainScriptPath <- "~/Desktop/TreeD"
setwd(mainScriptPath)
lapply(list.files("libs",full.names=T,recursive=T), source )
# Load RGL
require(rgl)



# Take in a phylogeny
treefile <- file.path(mainScriptPath, "chloro_anno2.tree")
phy <- readCommentedNexusTree(treefile[1]) #read 1st one only



#	If I get values for the tips in a matrix, we can estimate the ancestors via ML with ace() and declare that as my zz coords.
#	zz <- rnorm(Ntip(phy),1,0.15)
#	zz <- c(zz,ace(zz,phy)$ace)

# Alternativelly, I can get the information form the rate parameter from the tree for instance:

diffRate <- getDeltaMolEvoRate(phy,phy@data$rate)
accum <- getAccumSubsAtNode(phy,diffRate)

# I am not sure of how to deal with the ZZ coordinates now. They might be best stored as data in a phylo4d
# object. Need to check what option will hold well when reordering branches, ladderizong, etc:
# phy <- phylo4d(phy,all.data=zz)

# For now I'll just calculte the xy coords for the nodes and add ZZ coords to the obj.
xy <- xyPhyloCoord(phy)
zz <- accum
xyz <- cbind(xy,zz) # The ploting methods like the name zz. FIX THAT, so ANY vector name can be given to


# Ploting methods

# To do the actual plotting in RGL we:
# Open a new glwindow (in x11 for mac)
rgl::open3d()
rgl::rgl.material(line_antialias =F)
rgl::bg3d(col="black")

# Plot the nodes as spheres and edges as lines, so this will look like a "balls&sticks" style...
# Ploting branches as cylinders might be interesting too.
# For simplicity I am passing single values for attributes like shpere radius, line thickness and
# color, but these could be vectors too, so one can make line thickness ~ to support, pop size, etc.

mycolors <- c(rep("grey",phylobase::nTips(phy)),rep("dodgerblue4",phylobase::nNodes(phy)))
#mycolors = c(rep("aqua",Ntip(phy)),rep("black",Nnode(phy))) 				# APE VERSION

plotNodes3d(xyz,radius=.015,col= mycolors)
aspect3d(1.5,1.5,.8)


#plotEdges3d(phy,xyz,style="square", edgeType="cylinder",cylRadius=0.01,lineLwd=5,col="grey35")
plotEdges3d(phy,xyz,style="square", edgeType="line",cylRadius=0.005,lineLwd=5,col="dodgerblue4")


# Change the aspect (xyz proportions) of the plot. I can bind that to the parent window size or
#aspect3d(1.5,1.5,.8)

# Here we draw grids along on the Z and X end planes... Makes easier to visualize the info.
rgl::grid3d("z",col="white",lwd=0.5,n=5)
rgl::grid3d("x+",col=c("white","grey"),lwd=0.5,n=10)
#	rgl::axes3d("x+",cex=1,col=c("grey35"),tick = TRUE)
#	rgl::title3d(xlab="Time",zlab=NULL,col="white",cex=1.5)


# Labeling is an issue. For visualizing on the computer screen we want labels facing the user, that is, 
# facing the positive Z axis despite of the camera moviments. Loabeling a physical object would be tricky
# though, need to put more thought into it.
#	rgl::title3d(xlab="Time",zlab=NULL,col="white",cex=1.5)
#	plotLabels3d(phy,xyz,color="gray",cex=.7,adj=c(0,0))

rgl.postscript(filename="~/Blah.pdf",fmt="pdf")

# If you want to change the focus to a specific rgl device, run this cmd with the device's number
rgl::rgl.set(1)

##########################################################################################################
