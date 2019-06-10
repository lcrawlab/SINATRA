####### Generate Directions #######

#'Generate Equidistributed Cones
#'
#' @export
#'@description  Generate cones that are equidistributed about the sphere. We first generate equidistributed points, then use the rodrigues angle formula
#'to compute cones about those points.
#'
#'@param num_directions (int): The number of equidistributed directions we want on the sphere.
#'@param cap_radius (float): The radius of the cones we generate (determines the size of each cone).
#'@param directions_per_cone (int): The number of directions we want generated within each cone.
#'
#'@return directions (num_directions*directions_per_cone x 3 matrix): A matrix of equidistributed cones.
generate_equidistributed_cones <- function(num_directions, cap_radius, directions_per_cone){
  # generate central directions that are equidistributed around the sphere.
  cones <- generate_equidistributed_points(num_directions, num_directions)

  # renormalize these directions
  cones <- t(apply(cones, 1, function(x) x/sqrt(sum(x*x))))


  # generate directions for each cone
  directions <- matrix(0,ncol=3,nrow=0)
  for (i in 1:(dim(cones)[1]) ){
    directions <- rbind(directions,cones[i,])
    if (directions_per_cone > 1){
      directions <- rbind(directions,rodriq(cones[i,],cap_radius,directions_per_cone-1))
    }
  }

  directions
}
generate_random_cones <- function(num_directions,cap_radius,directions_per_cone){
  # uniformly generate random directions on the sphere
  cones <- matrix(rnorm(3*num_directions),ncol = 3)

  # renormalize these directions
  cones <- t(apply(cones, 1, function(x) x/sqrt(sum(x*x))))

  # generate directions for each cone
  directions <- matrix(0,ncol=3,nrow=0)
  for (i in 1:num_directions){
    directions <- rbind(directions,rodriq(cones[i,],cap_radius,directions_per_cone))
  }


  directions
}

#' Generate Equidistributed points on a sphere
#'
#' @export
#' @description Generate equidistributed points about a sphere using (insert reference)
#'
#' @param desired_number (int): the desired number of equidistributed points on the 2-sphere.
#' @param N (int): the initial number of points that the algorithm will try to generate. If the number of points generated is less than the desired number, the function will increment N
#' @return points (desired_number x 3 matrix): The matrix of equidistributed points in on the sphere.
generate_equidistributed_points <- function(desired_number, N){
  a <- 4*pi/N
  d <- sqrt(a)
  points <- matrix(NA,0,3)

  M_theta <- round(pi/d)
  d_theta <- pi/M_theta
  d_phi <- a/d_theta
  for (i in 0:(M_theta-1)){
    theta <- pi*(i + 0.5)/M_theta
    M_phi <- round(2*pi*sin(theta)/d_phi)
    for (j in 0:(M_phi - 1)){
      phi <- 2*pi*j/M_phi
      point <- c( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
      points <- rbind(points,point)
    }
  }
  if (dim(points)[1] < desired_number){
    return(generate_equidistributed_points(desired_number,N+1))
  }
  else{
    points = matrix(points[1:desired_number,],ncol = 3)
    return(points)
  }
}

#'
#' Compute the directions about a cone using the rodrigues angle formula.
#'  @export
#' @param z (vector of length 3) The direction around which we want to rotate.
#' @param r (float): The radius of the cone we want. Controls the cone size.
#' @param j (int): The number of directions in the cone.
#' @return dir (jx3 matrix): of directions in a cone of radius r about z.
rodriq<-function(z,r,j){
  z<-z/sqrt(z[1]^2+z[2]^2+z[3]^2)
  if( z[1]*z[2]*z[3]==0){
    z0<-c(0,0,0)
    z0<-c(1,1,1)*(z[]==0)
  }
  else{
    z0<-z
    z0[1]<- 2/z[1]
    z0[2]<- -1/z[2]
    z0[3]<- -1/z[3]
  }
  z0<-r*(z0/sqrt(z0[1]^2+z0[2]^2+z0[3]^2))
  z0<-z+z0
  z0<-z0/sqrt(z0[1]^2+z0[2]^2+z0[3]^2)
  B<-cross(z,z0)
  C<-z*(as.vector(z%*%z0))
  #print(B)
  #print(C)
  dir<-matrix(0,ncol=3,nrow=j)
  for (i in 1:j) {
    dir[i,]<-z0*cos(2*pi*i/j)+B*sin(2*pi*i/j)+C*(1-cos(2*pi*i/j))

  }
  return(dir)
}

#' Cross Product
#'
#' @description Compute the Cross Product
#' @param x (vector): the first (1,3) vector in the cross product.
#' @param y (vector): the second (1,3) vector in the cross product.
#' @return a (vector): the cross product of x and y.
cross<-function(x,y){
  a<-x
  a[1]<-x[2]*y[3]-x[3]*y[2]
  a[2]<-x[3]*y[1]-x[1]*y[3]
  a[3]<-x[1]*y[2]-x[2]*y[1]
  return(a)
}
