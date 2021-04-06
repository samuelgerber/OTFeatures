#' Title
#' Optimal Transport Allocation and Transport between 3D images
#'
#' Compute optimal transport between X and Y and create allocation
#' and transport images
#'
#'
#' @param from.image 3d array
#' @param to.image 3d array
#' @param mass.cost cost of allocating or removing mass
#' @param transport.type  optimal transport mode [default=3]: \n
#'    0 - Balanced (equalizes source and target mass)
#'    1 - Add mass (Only allow addition of mass in soutrce or target)
#'    2 - Subtract mass (Only allow subtraction of mass in source or target)
#'    3 - Match source mass (Add or subtract in source to match target)
#'    4 - Free (Add or subtract mass anywhere)
#'    5 - Free source only (Add or subtract mass in source only)
#' @param p.degere optimal transport cost degree (i.e. cost^p.degreee)
#'
#' @return list with allocation and trasnport cost images in from and to
#'
#' @examples
#' data(brains)
#' otf = extract.otf.image.3d(brain1, brain2)
#' par(mfrow=c(2,2))
#' image(brain1[,,2])
#' title("From")
#' image(brain2[,,2])
#' title("To")
#' image(brain1[,,2] - brain2[,,2])
#' title("Difference Per Voxel")
#' image(otf$difference.from[,,2])
#' title("Mass Allocation in From")
#'
#' @export extract.otf.image.3d
extract.otf.image.3d <- function(from.image, to.image, mass.cost=0, transport.type=3, p.degree=2) {
  XP = which(from.image>0, arr.ind=TRUE)
  XW = from.image[XP]
  YP = which(to.image>0, arr.ind=TRUE)
  YW = to.image[XP]
  trp = run.otf.points(XP, YP, XW, YW, mass.cost, transport.type, p.degree)

  k <- length(trp$cost)
  fromMass = trp$fromMass[[k]]
  toMass = trp$toMass[[k]]
  fromCoords = trp$from[[k]]
  toCoords = trp$to[[k]]
  map = trp$map[[k]]

  from <- as.vector(fromMass)
  to   <- as.vector(toMass)
  from.total <- sum( from )
  to.total   <- sum( to )

  library(data.table)
  map <- data.table::as.data.table(map)
  map.from <- map[,.( mass=sum(mass), cost = sum(cost*mass)/sum(mass) ), by=from]
  map.to <- map[,.( mass=sum(mass), cost = sum(cost*mass)/sum(mass) ), by=to]
  to[map.to[,to]] = to[map.to[,to]] - map.to[, mass]
  from[map.from[,from]] = from[map.from[,from]] - map.from[, mass]


  dim.x <- dim(from.image)[1]
  dim.y <- dim(from.image)[2]
  dim.z <- dim(from.image)[3]


  difference.from = array( 0, dim = c(dim.x, dim.y, dim.z) )
  transport.from = array( 0, dim = c(dim.x, dim.y, dim.z) )
  from <- data.table::data.table( x=fromCoords[,1], y=fromCoords[,2], z=fromCoords[,3] , from = from )
  from <- from[, .(from = sum(from)), by=.(x,y,z)]
  difference.from[ as.matrix(from[, .(x, y, z)]) ] = from[, from]

  from.cost <- numeric( length(fromMass) )
  from.cost[ map.from[,from]] <- map.from[,cost]
  from.cost <- data.table::data.table( x=fromCoords[,1], y=fromCoords[,2], z=fromCoords[,3], cost=from.cost)
  from.cost <- from.cost[, .(cost=sum(cost)), by=.(x,y,z)]
  transport.from[ as.matrix( from.cost[, .(x, y, z)]) ] = from.cost[, cost]

  dim.x <- dim(to.image)[1]
  dim.y <- dim(to.image)[2]
  dim.z <- dim(to.image)[3]

  difference.to = array( 0, dim = c(dim.x, dim.y, dim.z) )
  transport.to = array( 0, dim = c(dim.x, dim.y, dim.z) )

  to <- data.table::data.table( x=toCoords[,1], y=toCoords[,2], z=toCoords[,3], to = to )
  to <- to[, .(to=sum(to)), by=.(x,y,z)]
  difference.to[ as.matrix( to[, .(x, y, z)] ) ] = to[, to]

  to.cost <- numeric( length(toMass) )
  to.cost[map.to[,to]] <- map.to[,cost]
  to.cost <- data.table::data.table( x=toCoords[,1], y=toCoords[,2], z=toCoords[,3], cost=to.cost)
  to.cost <- to.cost[, .( cost = sum(cost)), by=.(x,y,z)]
  transport.to[ as.matrix(to.cost[, .(x, y, z)]) ] = to.cost[, cost]


  difference.to[ abs(difference.to) < 0.001 * min(c(toMass, fromMass)) ] = 0
  difference.from[ abs(difference.from) < 0.001 * min( c(toMass,fromMass) ) ] = 0
  #difference.to[ abs(difference.to) < 10^-10 ] = 0
  #difference.from[ abs(difference.from) < 10^-10 ] = 0

  #list( from.total=from.total, to.total = to.total, to=to, from=from,
  list( minimum = min(from.total, to.total), difference.to=difference.to,
        difference.from=difference.from,
        transport.from=transport.from, transport.to=transport.to)

}


#' Title
#' Optimal Transport Allocation and Transport between Points
#'
#' Compute optimal transport between X and Y and create allocation
#' and transport images
#'
#'
#' @param from 3d array
#' @param to 3d array
#' @param from.weights weights of from points
#' @param to.weights weights of to points
#' @param mass.cost cost of allocating or removing mass
#' @param transport.type  optimal transport mode [default=3]: \n
#'    0 - Balanced (equalizes source and target mass)
#'    1 - Add mass (Only allow addition of mass in soutrce or target)
#'    2 - Subtract mass (Only allow subtraction of mass in source or target)
#'    3 - Match source mass (Add or subtract in source to match target)
#'    4 - Free (Add or subtract mass anywhere)
#'    5 - Free source only (Add or subtract mass in source only)
#' @param p.degere optimal transport cost degree (i.e. cost^p.degreee)
#'
#' @return list woth allocation and trasnport cost images in from and to
#'
#' @export extract.otf.points
extract.otf.points <- function(from, to, from.weight=1, to.weight=1, mass.cost=0, transport.type=3, p.degree=2) {

  trp = run.otf.points(from, to, from.weight, to.weght, mass.cost, transport.type, p.degree)

  k <- length(trp$cost)
  fromMass = trp$fromMass[[k]]
  toMass = trp$toMass[[k]]
  fromCoords = trp$from[[k]]
  toCoords = trp$to[[k]]
  map = trp$map[[k]]

  from <- as.vector(fromMass)
  to   <- as.vector(toMass)
  from.total <- sum( from )
  to.total   <- sum( to )

  map <- data.table::as.data.table(map)
  map.from <- map[,.( mass=sum(mass), cost = sum(cost*mass)/sum(mass) ), by=from]
  map.to <- map[,.( mass=sum(mass), cost = sum(cost*mass)/sum(mass) ), by=to]
  to[map.to[,to]] = to[map.to[,to]] - map.to[, mass]
  from[map.from[,from]] = from[map.from[,from]] - map.from[, mass]


  from <- data.table::data.table( (coords):=fromCoords , from = from )
  from <- from[, .(from = sum(from)), by=.(coords)]

  from.cost <- numeric( length(fromMass) )
  from.cost[ map.from[,from]] <- map.from[,cost]
  from.cost <- data.table::data.table( (coords):=fromCoords , cost=from.cost)
  from.cost <- from.cost[, .(cost=sum(cost)), by=.(coords)]


  to <- data.table::data.table( (coords):=toCoords , to = to )
  to <- to[, .(to=sum(to)), by=.(coords)]

  to.cost <- numeric( length(toMass) )
  to.cost[map.to[,to]] <- map.to[,cost]
  to.cost <- data.table::data.table( (coords):=toCoords , cost=to.cost)
  to.cost <- to.cost[, .(cost=sum(cost)), by=.(coords)]


  list(from.mass=from, from.cost=from.cost, to.mass=to.mass, to.cost=to.cost)

}

run.otf.points <- function(X, Y, WX=1, WY=1, mass.cost=0, transport.type=2, p.degree=2) {

  trp.lp <- mop::multiscale.transport.create.lp( oType = 31, transport.type=transport.type,
                                          massCost=mass.cost )
  icprop <- mop::multiscale.transport.create.iterated.capacity.propagation.strategy( 3, 0 )
  mop::multiscale.transport.set.propagation.strategy.1( trp.lp, icprop )

  gmra1 <- gmra::gmra.create.ikm(X = X, eps = 0, nKids=8, stop=3)
  gmra2 <- gmra::gmra.create.ikm(X = Y, eps = 0, nKids=8, stop=3)
  trp <- mop::multiscale.transport.solve( trp.lp, gmra1, gmra2, p = p.degree, nType=1,
                                          dType=1, scaleMass = transport.type==0,
                                          w1=WX, w2 = WY )
  trp
}
