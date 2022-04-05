require(ggplot2)
## Loading required package: ggplot2
## Warning: package 'ggplot2' was built under R version 3.1.2
#set up the data frame.  The outer ring is defined in an anti-clockwise direction and the
#inner holes clockwise.  I have also labelled the holes.



holePolygon <- function(pos){
  
  
  pos[pos$hole == '0', 5] <- 'FALSE'
  pos[pos$hole == '1', 5] <- 'TRUE'
  
  
  
  
  closest_to_first_row <- function (xy){
    if (nrow(xy)==2){
      return (2)
    } else {
      return(1 + which.min(rowSums((xy[-1,] - xy[1,])^2)))
    }
  }
  
  bridge <- function(xy1, xy2){
    #which row in xy2 is closeset to mid point of xy1
    i.2 <- which.min(apply((colMeans(apply(xy1, 2,range))-xy2)^2,1,sum))
    i.1 <- which.min(apply((xy2[i.2,]-xy1)^2,1,sum)) #which row in xy1 is closest to this point in xy2
    i.2 <- which.min(apply((xy1[i.1,]-xy2)^2,1,sum)) #and again
    return(as.integer(c(i.1,i.2)))
  }
  
  require(data.table)
  ## Loading required package: data.table
  pos<-as.data.table(pos)
  
  #Add a new field to contain the original position as I need this for
  #drawing the orginal borders (which are correct)
  pos$draw=pos$ring
  
  while(TRUE){
    #summarise the information on each ring
    holes <- pos[hole==TRUE,.(mid.x=mean(range(x)), mid.y=mean(range(y))),by=.(id,ring)]
    holes <- holes[,num_holes:=.N ,by=id]
    if(max(holes[,num_holes])<=1) break  #Exit if only one hole per id remaining
    
    #get the first id that has more than one hole in it
    ii <- which(holes$num_holes>1)[1]
    h.id <- holes$id[ii]                  #which id are we dealing with           
    h1.ring <- holes$ring[ii]             #which hole are we dealing with first
    
    #get the hole ring which is closest to h1.ring  This ensures that the shortest path 
    #between h1.ring and h2.ring does not cross another hole ring.
    h2.ring <- holes[id==h.id][closest_to_first_row(as.matrix(holes[id==h.id, .(mid.x,mid.y)])),ring]
    cat('joining ring', h1.ring, 'to ring', h2.ring, '\n')
    
    #find the best bridging point
    h1.xy <- as.matrix(pos[id==h.id & ring==h1.ring, .(x, y)])             #xy matrix for ring1
    h2.xy <- as.matrix(pos[id==h.id & ring==h2.ring, .(x, y)])             #xy matrix for ring2
    h1.l  <- nrow(h1.xy)                                                   #number of points in ring1
    h2.l  <- nrow(h2.xy)                                                   #number of points in ring2
    h1.draw <- pos[id==h.id & ring==h1.ring, draw]                         #existing values for drawing border
    h2.draw <- pos[id==h.id & ring==h2.ring, draw]
    
    b <- bridge(h1.xy, h2.xy)   #b[1] is the row in h1 and b[2] is the row in h2 to bridge
    
    #reorder h2 values about the bridging point and insert into the bridge point in h1
    new.xy <-   rbind(
      h1.xy[seq(b[1]),]            #h1 points up to the bridge
      ,h2.xy[seq(b[2], h2.l-1),]    #h2 from over the bridge to one before the tail=head
      ,h2.xy[seq(1,b[2]),]          #h2 from the head to the bridge again
      ,h1.xy[seq(b[1], h1.l),]      #h1 from the bridge to the tail
    )
    new.draw <- c( h1.draw[seq(b[1])]            #arrange the 'draw' to line up with the orginal rings
                   ,h2.draw[seq(b[2], h2.l-1)]   #so can jump from one ring to another without drawing
                   ,h2.draw[seq(1,b[2])]         #a border over the jump
                   ,h1.draw[seq(b[1], h1.l)]  
    )
    
    #delete the old values and replace with the new values
    drop.rows <- which(pos$id==h.id & (pos$ring==h1.ring|pos$ring==h2.ring))
    
    #update the pos data frame by dropping the original and adding the new
    pos <- rbind(pos[-drop.rows,]
                 ,data.frame(id=h.id
                             ,ring=h1.ring
                             ,x=new.xy[,1]
                             ,y=new.xy[,2]
                             ,hole=TRUE
                             ,draw=new.draw)
    )
  }
  
  #reorder the pos data.frame according to the new rings (with the holes merged)
  pos<-pos[order(id,ring),]
  
  
  # return
  bdry_tr <- pos
}

Network <- function(cellPos, prox_thresh){
  
  r <- tri.mesh(cellPos$X, cellPos$Y)
  delaunayTri <- tripack::triangles(r)
  
  Num_tri <- nrow(delaunayTri)
  # coord length of triangle
  a <- delaunayTri[,1] # node 1 index
  b <- delaunayTri[,2] # node 2 index
  c <- delaunayTri[,3] # node 3 index
  tri_sidelength <- cbind(sqrt( ( r$x[a] - r$x[b] ) ^ 2 + ( r$y[a] - r$y[b] ) ^ 2 ), sqrt( ( r$x[a] - r$x[c] ) ^ 2 + ( r$y[a] - r$y[c] ) ^ 2 ), sqrt( ( r$x[b] - r$x[c] ) ^ 2 + ( r$y[b] - r$y[c] ) ^ 2 ))
  
  
  # perimeter of triangles
  tri_perimeter <- rowSums(tri_sidelength)
  
  # area of triangles
  s <- 0.5 * tri_perimeter
  tri_area <- sqrt(s * (s - tri_sidelength[, 1]) * (s - tri_sidelength[, 2]) *
                     (s - tri_sidelength[, 3]))
  
  # clustering
  k <- seq_len( r$tlnew - 1 )
  i <- r$tlist[k]
  j <- r$tlist[r$tlptr[k]]
  keep <- i > 0
  
  
  i <- abs(i[keep])
  j <- abs(j[keep])
  
  
  distances <- sqrt( ( r$x[i] - r$x[j] ) ^ 2 + ( r$y[i] - r$y[j] ) ^ 2 )
  
  i <- i[distances <= prox_thresh]
  j <- j[distances <= prox_thresh]
  #--------------------------------------------#
  
  Dual_NodeList <- cbind(seq(1, nrow(cellPos)),  cellPos[,c('X', 'Y')])
  colnames(Dual_NodeList) <- c('nodes', 'x', 'y')
  
  
  Dual_EdgeList <- data.frame(cbind(i, j))
  colnames(Dual_EdgeList) <- c('col1', 'col2')
  
  Dual_EdgeList <- Dual_EdgeList %>%
    mutate(from = pmin(col1, col2),
           to = pmax(col1, col2)) %>%
    distinct(from, to)
  
  
  colnames(Dual_EdgeList) <- c('from', 'to')
  
  
  # return lists
  return(list(Dual_EdgeList, Dual_NodeList))
  
}



poisp <- function(pos_Target, rect, nsim){
  
  #pos_Target <- roiDat[roiDat$class == typegrid[id,1], ]
  #nsim <- 500
  n <- nrow(pos_Target)
  
  if(missing(rect)){
    xmax <- max(pos_Target$X)
    ymax <- max(pos_Target$Y)
    
    xmin <- min(pos_Target$X)
    ymin <- min(pos_Target$Y)
  } else {
    xmax <- max(rect$X)
    ymax <- max(rect$Y)
    
    xmin <- min(rect$X)
    ymin <- min(rect$Y)
  }
  #plot(NeighborDat)
  simpp <- runifpoint(n = n, win = owin(c(xmin, xmax), c(ymin, ymax)), nsim = nsim)
  #plot(simpp$`Simulation 2`$x, simpp$`Simulation 2`$y)
  # return the result 
  return(simpp)
}


#-------------------------------------------------------#
# Find number of interactions, voronoi neighbor verison #
#-------------------------------------------------------#


nnIntrxn.voi <- function(roiDat, type1, type2, dThresh, nsim){
  
  #type1 <- ctype1
  #type2 <- ctype2
  #dThresh <- 75
  #nsim <- 500
  # margin
  Rect <- data.frame(cbind(X = c(min(roiDat$X), min(roiDat$X), max(roiDat$X), max(roiDat$X)), Y = c(min(roiDat$Y), max(roiDat$Y), max(roiDat$Y), min(roiDat$Y))))
  
  # randomize - 1
  #set.seed(999)
  simpp_c1 <- poisp(roiDat[roiDat$class == as.character(type1), ], rect = Rect, nsim = nsim) # randomize cell type 1's location
  
  #set.seed(999)
  simpp_c2 <- poisp(roiDat[roiDat$class == as.character(type2), ], rect = Rect, nsim = nsim) # randomize cell type 2's location
  
  # loop through list
  nulldist <- foreach(i = 1:500, .combine = 'c', .packages = c('RANN', 'reshape2')) %dopar% {
    
    # scenairo 1: randomized 1 versus true 2
    cdata1 <- cbind(simpp_c1[[i]]$x, simpp_c1[[i]]$y)
    
    
    dists <- nn2(cdata1, query = roiDat[roiDat$class == as.character(type2), c('X', 'Y')], k = nrow(cdata1), treetype = 'kd', searchtype = 'radius', radius = 78)
    #flexclust::dist2(cdata1, roiDat[roiDat$class == as.character(type2), c('X', 'Y')]) # cdata1 - row
    
    nn.idx <- data.frame(dists$nn.idx)
    
    nn.idx$source <- seq_len(nrow(roiDat[roiDat$class == as.character(type2), c('X', 'Y')]))
    
    nn.merge <- melt(nn.idx, id.vars = 'source')
    
    # clear non-valid rows
    nn.merge <- nn.merge[nn.merge$value != 0, -2]
    
    
    simCount1 <- nrow(nn.merge)
    #df1 <- cbind.data.frame(seq_len(nrow(cdata1)), apply(dists, 1, min)) %>%
    #  `colnames<-`(c("row.id", 'min.dist')) %>%
    #  filter(min.dist <= dThresh) %>%
    #  nrow()
    
    
    # scenairo 2: randomized 2 versus true 1
    cdata2 <- cbind(simpp_c2[[i]]$x, simpp_c2[[i]]$y)
    #dists <- flexclust::dist2(roiDat[roiDat$class == as.character(type1), c('X', 'Y')], cdata2) # cdata1 - row
    #df2 <- cbind.data.frame(seq_len(nrow(roiDat[roiDat$class == as.character(type1), c('X', 'Y')])), apply(dists, 1, min)) %>%
    #  `colnames<-`(c("row.id", 'min.dist')) %>%
    #  filter(min.dist <= dThresh) %>%
    #  nrow()
    dists <- nn2(cdata2, query = roiDat[roiDat$class == as.character(type1), c('X', 'Y')], k = nrow(cdata2), treetype = 'kd', searchtype = 'radius', radius = 78)
    #flexclust::dist2(cdata1, roiDat[roiDat$class == as.character(type2), c('X', 'Y')]) # cdata1 - row
    
    nn.idx <- data.frame(dists$nn.idx)
    
    nn.idx$source <- seq_len(nrow(roiDat[roiDat$class == as.character(type1), c('X', 'Y')]))
    
    nn.merge <- melt(nn.idx, id.vars = 'source')
    
    # clear non-valid rows
    nn.merge <- nn.merge[nn.merge$value != 0, -2]
    
    simCount2 <- nrow(nn.merge)
    
    #nulldist <- c(nulldist, 0.5*(df1 + df2))
    simnull <- 0.5 * (simCount1 + simCount2)
  }
  return(nulldist)
}



#-----------------------------#
# Find number of interactions #
#-----------------------------#

nnIntrxn <- function(sDF, simDF, radius){ # sDF: the coordiantes for the source point pattern
  
  #sDF <- sPos
  if(nrow(sDF) == 0 | nrow(simDF) == 0){
    simInt <- 0
  } else{
    dist <- nn2(simDF, query = sDF, k = nrow(simDF), treetype = 'kd', searchtype = 'radius', radius = radius)
    
    nn.idx <- data.frame(dist$nn.idx)
    
    nn.idx$source <- seq_len(nrow(sDF))
    
    nn.merge <- reshape2::melt(nn.idx, id.vars = 'source')
    
    # clear non-valid rows
    nn.merge <- nn.merge[nn.merge$value != 0, -2]
    
    
    # remove itself
    
    
    nn.merge$diff <- nn.merge$source - nn.merge$value
    nn.merge <- nn.merge[nn.merge$diff != 0, -3]
    
    simInt <- nrow(nn.merge)
  }
  
  
  # replace by real cell type
  return(simInt)
}

#---------------------------------------#
#-------- Bivariate Kcross Function ----#
#---------------------------------------#


bivarAnalysis.Kcross <- function(pts.type1, pts.type2, Region){
  
  #Region <- Rect
  #pts.type1 <- 'posCD4'
  #pts.type2 <- 'posCD8'
  
  type1 <- get(eval(pts.type1)) %>%
    select('X', 'Y')
  
  type2 <- get(eval(pts.type2)) %>%
    select('X', 'Y')
  
  if(nrow(type1)*nrow(type2) != 0){
    
    # read pts dat
    #Region <- Region_HE
    
    type1$attr <- pts.type1
    
    type2$attr <- pts.type2
    
    # create multitype df
    pts_OI <- rbind(type1, type2)
    
    # define the type
    species <- factor(pts_OI$attr)
    
    # create multitype ppp
    #Region <- Region_CK56
    
    
    # check if empty  
    ppp1 <- ppp(type1$X, type1$Y, owin(poly = Region))
    ppp2 <- ppp(type2$X, type2$Y, owin(poly = Region))
    
    # prevent NA 
    
    if(is.empty(ppp1) == 'FALSE' & is.empty(ppp2) == 'FALSE'){
      
      
      
      
      multitype_ppp <- ppp(pts_OI$X, pts_OI$Y, marks = species, owin(poly = Region))
      K.cross <- data.frame(Gcross(multitype_ppp, i = pts.type1, j = pts.type2, r = seq(0,78,0.1), correction = 'border'))
      
      #plot(Gihc)
      # relocat DF
      
      K.cross <- K.cross[complete.cases(K.cross),]
      
      
      
      # calculate the area (positive - negative )  
      K.cross$km <- K.cross$km - K.cross$theo
      
      i.to.j.diff.area <- trapz(K.cross$r, K.cross$km) 
      
      
      # j to i
      
      K.cross <- data.frame(Gcross(multitype_ppp, i = pts.type2, j = pts.type1, r = seq(0,78,0.1), , correction = 'border'))
      
      K.cross <- K.cross[complete.cases(K.cross),]
      
      
      
      K.cross$km <- K.cross$km - K.cross$theo
      K.cross <- K.cross[complete.cases(K.cross),]
      j.to.i.diff.area <- trapz(K.cross$r, K.cross$km) 
      
      
    } else {
      
      i.to.j.diff.area <- 0
      j.to.i.diff.area <- 0
      
    }
  } else {
    i.to.j.diff.area <- 0
    j.to.i.diff.area <- 0
  }
  return(list(i.to.j.diff.area, j.to.i.diff.area))
}






#------- Shannon Entropy --------------#

ShannonE <- function(types, coreDat){
  
  #types <- ctype_tensor
  #coreDat <- ROIDat
  type.count <- length(types)
  
  Total <- nrow(coreDat) # total number of cells
  
  # all int/ext data
  ctype_stat_all <- data.frame(matrix(nrow =  0, ncol = 0))
  
  for (cseq in seq_len(type.count)) {   # outer loop, calcualte interior stats
    #cseq <- 1
    ctype_int <- 0
    ctype_ext <- 0
    p <- 0 # ratio
    
    # current cell type
    ctype <- types[cseq]
    
    # coordinates data for the current core, current cell type
    
    ctype.Dat <- coreDat[coreDat$class == ctype, c('X', 'Y')]
    
    p <- nrow(ctype.Dat)/Total
    
    # interior score
    ctype_int <- mean(as.matrix(dist(ctype.Dat[,c('X', 'Y')])))
    
    
    # other cell types
    ctype_other <- types[-cseq]
    
    ctype_stat <- cbind(p, ctype_int)
    
    for (cseq_other in ctype_other) {
      
      ctype_other.Dat <- coreDat[coreDat$class == cseq_other, c('X', 'Y')]
      
      # exterior score
      
      ctype_ext <- ctype_ext +  mean(flexclust::dist2(ctype.Dat, ctype_other.Dat))
      
    }
    # number of computations = No. cell types - 1
    ctype_stat <- cbind(ctype_stat, ctype_ext / (type.count - 1))
    
    colnames(ctype_stat) <- c('p', 'int', 'ext')
    
    ctype_stat_all <- rbind(ctype_stat_all ,cbind(ctype_stat, ctype))
    
  }
  
  ShannonH <- 0
  # combine row data
  
  
  for (dat in seq_len(type.count)) {
    p <- as.numeric(as.character(ctype_stat_all[dat, 1]))
    d_int <- as.numeric(as.character(ctype_stat_all[dat, 2]))
    
    d_ext <- as.numeric(as.character(ctype_stat_all[dat, 3]))
    d_final <- d_int / d_ext
    
    if(isTRUE(d_int*d_ext == 0)){
      d_final <- 0
    }
    if(isTRUE(p != 0)){
      ShannonH <- -d_final*p*log2(p) + ShannonH
    }
  }
  
  return(ShannonH)
}




#---------------------------------------#
#-------- Bivariate Kcross Function ----#
#---------------------------------------#

gcross <- function(pts.type1, pts.type2, Region){
  
  #Region <- Rect
  #pts.type1 <- 'posCD4'
  #pts.type2 <- 'posCD8'
  
  type1 <- get(eval(pts.type1)) %>%
    select('X', 'Y')
  
  type2 <- get(eval(pts.type2)) %>%
    select('X', 'Y')
  
  
  # read pts dat
  #Region <- Region_HE
  
  type1$attr <- pts.type1
  
  type2$attr <- pts.type2
  
  # create multitype df
  pts_OI <- rbind(type1, type2)
  
  # define the type
  species <- factor(pts_OI$attr)
  
  # create multitype ppp
  
  # check if empty  
  ppp1 <- ppp(type1$X, type1$Y, owin(poly = Region))
  ppp2 <- ppp(type2$X, type2$Y, owin(poly = Region))
  
  
  multitype_ppp <- ppp(pts_OI$X, pts_OI$Y, marks = species, owin(poly = Region))
  G.cross1 <- data.frame(Gcross(multitype_ppp, i = pts.type1, j = pts.type2, r = seq(0,100,5), correction = 'border')) %>%
    select(r, rs)
  
  G.cross2 <- data.frame(Gcross(multitype_ppp, i = pts.type2, j = pts.type1, r = seq(0,100,5), correction = 'border')) %>%
    select(r, rs)
  #plot(Gihc)
  
  return(list(G.cross1, G.cross2))
}




#+++++++++++++++++++++++++++++++++++++#
#++++++++ Enrichment score +++++++++++#
#+++++++++++++++++++++++++++++++++++++#

ratio.sim <- function(CD4T, Mpg, Rect, nsim){
  
  set.seed(999)
  #nsim <- 500
  IRPvars <- c("PD1.", "PDL1.", 'KI67.', 'EOMES.',  "IL10.", 'ICOS.')
  
  
  simpp_CD4 <- poisp(CD4T, Rect, nsim = nsim) # randomize cell type 1's location
  
  simpp_Mpg <- poisp(Mpg, Rect, nsim = nsim) # randomize cell type 2's location
  
  ratios_CD4T <- data.frame(matrix(nrow = 0, ncol = 0 ))
  ratios_Mpg <- data.frame(matrix(nrow = 0, ncol = 0 ))
  
  for(sim in seq_len(500)){
    
    # current simulated CD4
    sim.CD4 <- cbind.data.frame(simpp_CD4[[sim]]$x, simpp_CD4[[sim]]$y) %>%
      `colnames<-` (c('X', 'Y')) %>%
      cbind.data.frame(CD4T[, 30:36])
    
    sim.Mpg <- cbind.data.frame(simpp_Mpg[[sim]]$x, simpp_Mpg[[sim]]$y) %>%
      `colnames<-` (c('X', 'Y')) %>%
      cbind.data.frame(Mpg[, 30:36])
    
    
    
    #set.seed(999)
    
    ### Case 1: CD4 T cell to Macrophage
    
    distMatrix <- flexclust::dist2(sim.Mpg[, c('X', 'Y')], CD4T[, c('X', 'Y')]) %>%
      as.matrix()
    
    # These CD4 T cells has at least 1 adjacent Macrophages
    adj.Mpg.id <- which(rowMins(distMatrix) < 100)
    
    # Get the single-cell data for these CD4 T cells
    adj.Mpg <- sim.Mpg[adj.Mpg.id, ]
    
    
    if(!(is.empty(adj.Mpg.id))){
      
      ratios <- adj.Mpg %>%
        summarise_each_(funs(sum(.==1,na.rm=TRUE)), IRPvars) %>%
        `/` (nrow(adj.Mpg))
      
      ratios$nsim <- sim
      
      ratios_Mpg <- rbind.data.frame(ratios_Mpg, ratios)
    }
    
    
    ### Case 2: Macrophage to CD4T
    
    distMatrix <- flexclust::dist2(sim.CD4[, c('X', 'Y')], Mpg[, c('X', 'Y')]) %>%
      as.matrix()
    
    # These CD4 T cells has at least 1 adjacent Macrophages
    adj.CD4T.id <- which(rowMins(distMatrix) < 100)
    
    # Get the single-cell data for these CD4 T cells
    adj.CD4T <- sim.CD4[adj.CD4T.id, ]
    
    
    if(!(is.empty(adj.CD4T.id))){
      
      ratios <- adj.CD4T %>%
        summarise_each_(funs(sum(.==1,na.rm=TRUE)), IRPvars) %>%
        `/` (nrow(adj.CD4T))
      
      ratios$nsim <- sim
      
      ratios_CD4T <- rbind.data.frame(ratios_CD4T, ratios)
    }
    
  }
  return(list(ratios_Mpg, ratios_CD4T))
  
}



#+++++++++++++++++++++++++++++++++++++#
#++++++++ Enrichment score - v2 +++++++++++#
#+++++++++++++++++++++++++++++++++++++#

ratio.sim2 <- function(CD4T, Mpg, Rect, nsim){
  
  set.seed(999)
  #nsim <- 500
  IRPvars <- c("PD1.", "PDL1.", 'KI67.', 'EOMES.',  "IL10.", 'ICOS.')
  
  
  simpp_CD4 <- poisp(CD4T, Rect, nsim = nsim) # randomize cell type 1's location
  
  simpp_Mpg <- poisp(Mpg, Rect, nsim = nsim) # randomize cell type 2's location
  
  ratios_CD4T <- data.frame(matrix(nrow = 0, ncol = 0 ))
  ratios_Mpg <- data.frame(matrix(nrow = 0, ncol = 0 ))
  
  for(sim in seq_len(500)){
    #sim <- 1  
    # current simulated CD4
    sim.CD4 <- cbind.data.frame(simpp_CD4[[sim]]$x, simpp_CD4[[sim]]$y) %>%
      `colnames<-` (c('X', 'Y')) %>%
      cbind.data.frame(CD4T[, c(30:36, 39)])
    
    sim.Mpg <- cbind.data.frame(simpp_Mpg[[sim]]$x, simpp_Mpg[[sim]]$y) %>%
      `colnames<-` (c('X', 'Y')) %>%
      cbind.data.frame(Mpg[, c(30:36, 39)])
    
    
    
    #set.seed(999)
    
    ### Case 1: Macrophage to CD4T
    
    distMatrix <- flexclust::dist2(sim.Mpg[, c('X', 'Y')], CD4T[, c('X', 'Y')]) %>%
      as.matrix()
    
    # These CD4 T cells has at least 1 adjacent Macrophages
    adj.Mpg.id <- which(rowMins(distMatrix) < 100)
    
    # Get the single-cell data for these CD4 T cells
    adj.Mpg <- sim.Mpg[adj.Mpg.id, ]
    
    
    if(!(is.empty(adj.Mpg.id))){
      
      ratios <- adj.Mpg %>%
        group_by(Status) %>%
        tally() %>%
        tidyr::complete(Status, fill = list(n = 0)) %>%
        mutate(ratio = n / nrow(adj.Mpg))  
      
      ratios$nsim <- sim
      
      ratios_Mpg <- rbind.data.frame(ratios_Mpg, ratios)
    }
    
    
    ### Case 2: CD4T to macrophage
    
    distMatrix <- flexclust::dist2(sim.CD4[, c('X', 'Y')], Mpg[, c('X', 'Y')]) %>%
      as.matrix()
    
    # These CD4 T cells has at least 1 adjacent Macrophages
    adj.CD4.id <- which(rowMins(distMatrix) < 100)
    
    # Get the single-cell data for these CD4 T cells
    adj.CD4 <- sim.CD4[adj.CD4.id, ]
    
    
    if(!(is.empty(adj.CD4.id))){
      
      ratios <- adj.CD4 %>%
        group_by(Status) %>%
        tally() %>%
        tidyr::complete(Status, fill = list(n = 0)) %>%
        mutate(ratio = n / nrow(adj.CD4))  
      
      ratios$nsim <- sim
      
      ratios_CD4T <- rbind.data.frame(ratios_CD4T, ratios)
    }
    
  }
  return(list(ratios_CD4T, ratios_Mpg))
  
}


setClass(Class="AshapePol",
         representation(
           goodAshape="logical",
           vert="numeric"
         )
)

getAlphaShapeInOrder <- function(ashapeX){
  returnVal=new("AshapePol", goodAshape=TRUE, vert=0)
  
  if (nrow(ashapeX$edges)==0) {
    #stop("Graph not connected")
    returnVal@goodAshape = FALSE
  }
  
  ashapeGraph = igraph::graph.edgelist(cbind(as.character(ashapeX$edges[, "ind1"]), 
                                             as.character(ashapeX$edges[, "ind2"])), directed = FALSE)
  
  
  
  #plot(ashapeGraph)
  ## check if is closed
  if (!igraph::is.connected(ashapeGraph)) {
    #stop("Graph not connected")
    returnVal@goodAshape = FALSE
  }
  if (any(igraph::degree(ashapeGraph) != 2)) {
    #stop("Graph not circular")
    returnVal@goodAshape = FALSE
  }
  if (igraph::clusters(ashapeGraph)$no > 1) {
    #stop("Graph composed of more than one circle")
    returnVal@goodAshape = FALSE
  }
  
  if(! returnVal@goodAshape){
    return(returnVal)
  }
  
  cutg = ashapeGraph - E(ashapeGraph)[1]
  # find chain end points
  ends = names(which(igraph::degree(cutg) == 1))
  path = get.shortest.paths(cutg, ends[1], ends[2])$vpath[[1]]
  # this is an index into the points
  pathX = as.numeric(V(ashapeGraph)[path]$name)
  # join the ends
  pathX = c(pathX, pathX[1])
  
  aShapePath <- Polygon(ashapeX$x[pathX, ])  
  # polygon(X[pathX, 1], X[pathX,2], col = "gray", border = "red")
  
  ## test within polygon
  inAshape = point.in.polygon(ashapeX$x[,1],ashapeX$x[,2],ashapeX$x[pathX,1],ashapeX$x[pathX,2])
  if (any(inAshape==0)){
    #stop("point outside alpha shape")
    returnVal@goodAshape = FALSE
  }else{
    returnVal@vert = pathX
  }
  return(returnVal)
}


getEllipse <- function(cls.data){
  
  alpha0 <- 1
  alphaStep <- 1

  goodAlpha = FALSE
  alpha <- alpha0
  #X<-subPoints[cidMatch,]
  X <- cls.data[, c('X', 'Y')]
  X <- unique(X)
  # try different alpha until shape covers all points
  while(!goodAlpha){
    print(alpha)
    alpha <- alpha + alphaStep
    ashapeX <- ashape(X, alpha = alpha)
    #print(plot(ashapeX, xlab = "x", ylab = "y"))

    res=getAlphaShapeInOrder(ashapeX)
    goodAlpha <- res@goodAshape
  }
  
  nPt <-nrow(X)
  
  ### get alpha-perimeter
  Pconc <- ashapeX$length
  
  ### get alpha-area:
  alphaShape <- Polygon(X[res@vert,])
  Aconc <- alphaShape@area

  hullPts <- chull(X)
  hullPts <- c(hullPts, hullPts[1])
  
  polygon_set <- matrix(nrow=nrow(X[hullPts,]),ncol=3)
  
  ### get convex hull
  
  #lines(X[hullPts, ])
  convShape <- Polygon(X[hullPts, ])  
  ## plot the polygon of convex hull !!!!!!!!
  
  polygon_set[,1] <- X[hullPts, 1]
  polygon_set[,2] <- X[hullPts, 2]
  polygon_set <- data.frame(polygon_set)
  
  ### convex hull area
  Aconv<-convShape@area
  
  ## fitting ellipse
  # we can get center and covariance matrix from the cluster points:
  ellFit <- cov.wt(X)
  # eigen vectors indicate orientations
  ellEigen <- eigen(ellFit$cov)
  # eigen values are associated with major/minor axes length (half)
  
  axes <- sqrt(ellEigen$values * qchisq(.95, df=2))
  
  # eccentricity
  eccentricity <- sqrt(1-(axes[2]/axes[1])^2) # 1: major, 2: minor
  
  # circularity
  circularity <- Aconc*4*pi/(Pconc)^2
  
  # convexity
  convexity <- Aconc/Aconv
  
  # assign label
  label <- 'else'
  if(eccentricity < 0.8 & convexity > 0.8 & circularity > 0.5){
    label <- 'circular'
  } 
  if(eccentricity > 0.9 | convexity < 0.3 | circularity < 0.3){
    label <- 'elongated'
  }
  
  # cluster morphometrics summary
  clus.morph <- cbind.data.frame(eccentricity, circularity, convexity, Aconv/1000000, label) %>%
    `colnames<-` (c('eccentricity', 'circularity', 'convexity', 'Aconv', 'label'))
  
}

