# for every direction, compute the euler chacteristic of a 3D complex, filtered by direction. 
#write dumb version that recomputes ec for every complex in the filtration
# more advanced version that computes the ec based on updates to the complex in the filtration.

# Complex (for simplicial complex) should contain vertices, faces, edges, and t? as in the EC3D.R file.
# We assume that Complex has a field Vertices such that complex$Vertices is a n x 3 matrix such that every row 
# gives the coordinates of the vertices. Complex$Edges should be a list of the edges, m x 2, where each row contains
# the vertices belonging to each edge. Complex$Faces has the same format.

# How to specify direction on the sphere? potentially use spherical coordinates, or unit vector coordinates.
# function should just be dot product of each point's coordinates with this unit vector, perhaps upon centering first.

# integrate the euler characteristic to get a curve.

library(geometry)
library(misc3d)
library(rgl)
library(Rvcg)
library(TDA)

# For point cloud data; one option for manual nearest neighbor type thing, or just use vietoris rips
# check https://stackoverflow.com/questions/22630226/3d-surface-plot-with-xyz-coordinates
create_simplicial_complex <- function(){
  # rely on alpha shape package; many other functions for surface reconstruction.
  # might need to play around with alpha
  ashape3d.obj <- ashape3d(x,alpha = 0.5)
  plot(ashape3d.obj) # to visualize
  
  vertices <- 1:length(x)
  edges <- ashape3d.obj$edge[which(ashape3d.obj$edge[,4] == 1),1:2] # column 4 of this object returns which edges are attached
  faces <- ashape3d.obj$triang[which(ashape3d.obj$triang[,4] == 1),1:3] 
  
  
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
}


# We assume the input is a text file in a similar format as an OFF file (not the same)
# The first line should say #vertices #faces #edges; the following lines are the vertices coordinates in each new line;
# following are the edges, with the indices of the vertices in the edges. The Following is are the vertices in each face, 
#preceded by the number of vertices in that face. See the example text file.
process_text_file <- function(fileName){
  con = file(fileName, "r") # con for connection
  lines = readLines(con) # apparently this is dangerous with big files
  close(con)
  split_line = strsplit(lines[1]," ")[[1]] #this is a list, we index the first element (a vector) with the double bracket 
  num_vertices = strtoi(split_line[1])
  num_edges = strtoi(split_line[3])
  num_faces = strtoi(split_line[2])
  n = num_vertices + num_edges + num_faces
  dim = length(strsplit(lines[2]," ")[[1]]) # get the dimension of the points
  face_size = strtoi(strsplit(lines[2 + num_vertices + num_edges], " ")[[1]])[1]
  
  
  # store the vertices
  vertices <- matrix(NA,nrow = num_vertices,ncol = dim)
  for (i in 2:(num_vertices+1)){
    vertices[i-1,] <- strtoi(strsplit(lines[i], " ")[[1]]) # coordinates of each vertex
  }
  
  #store the edges
  edges <- matrix(NA,nrow = num_edges,ncol = 2)
  for (i in (num_vertices+2):(num_vertices+num_edges+1)){
    edges[i - num_vertices - 1,] <- strtoi(strsplit(lines[i], " ")[[1]]) + 1# + 1 if the vertices are zero indexed
  }
  
  #store the faces
  faces <- matrix(NA, nrow = num_faces, ncol = face_size)
  for (i in (num_vertices+num_edges+2):(n+1)){
    faces[i - num_vertices - num_edges - 1,] <- 
      strtoi(strsplit(lines[i], " ")[[1]])[1:face_size+1]  +1 # + 1 if the vertices are zero indexed
  }
  
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
}
#Processing OFF files specifically

process_off_file_v2 <- function(fileName){
    con = file(fileName, "r") # con for connection
    lines = readLines(con) # apparently this is dangerous with big files
    close(con)
    split_line = strsplit(lines[2]," ")[[1]] #this is a list, we index the first element (a vector) with the double bracket 
    num_vertices = strtoi(split_line[1])
    num_edges = strtoi(split_line[3])
    num_faces = strtoi(split_line[2])
    n = num_vertices + num_faces
    dim = length(strsplit(lines[3]," ")[[1]]) # get the dimension of the points
    face_size = strtoi(strsplit(lines[3 + num_vertices], " ")[[1]])[1] # assumes faces have the same size
    
    
    # store the vertices and its coordinates
    vertices <- matrix(NA,nrow = num_vertices,ncol = dim)
    for (i in 3:(num_vertices+2)){
        vertices[i-2,] <- as.numeric(strsplit(lines[i], " ")[[1]]) # coordinates of each vertex
    }
    
    
    # edges to be read off below
    edges <- matrix(NA,nrow = 0,ncol = 2)
    
    #store the faces, process the edges at the same time
    faces <- matrix(NA, nrow = num_faces, ncol = face_size)
    for (i in (num_vertices+3):(n+2)){
        faces[i - num_vertices - 2,] <- 
            strtoi(strsplit(lines[i], " ")[[1]])[1:face_size+1] + 1 # + 1 if the vertices are zero indexed
        
        face_length = length(faces[i - num_vertices - 2,])
        # store edges in each face in a set; requires that faces list the vertices in increasing order.
        for (j in 1:(face_length-1) ){
            for (k in (j+1):(face_length) ){
                face_vertices <- sort(faces[i - num_vertices - 2,])
                edges <- rbind(edges, c(face_vertices[j],face_vertices[k]))
            }
        }
        
        
    }
    
    edges <- unique(edges)
    
    complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
}

process_off_file_v3=function(input_dir){
  off=vcgImport(input_dir,silent = TRUE)
  vertices=as.matrix(t(off$vb)[,1:3])
  faces=as.matrix(t(off$it))
  edges=vcgGetEdge(off)
  edges=as.matrix(edges[,1:2])
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
  return(complex)
}

# The complex is filtered by some function on the vertices of the simplicial complex; in this case, we use the distance
# to some plane in a certain direction in S2 (the two dimensional sphere). To compute this distance in R2, we simply multiply 
# coordinates of the vertices by a rotation matrix. In three dimensions, we dot product with the unit vector in the desired 
# direction.

#step size is how fine we want the ec curve to be.
# at the last step, we center the curve and integrate it to obtain the smooth euler characetristic curve transform (SECT)
# gEuler achieves this.

#The Euler characteristic of a subcomplex which in our case is a a surface of a polyhedra has a simple form: V - E + F
#curve_length outputs the ec curve as a vector of the same size. 
#vertex_function are the function values on the vertices
#SECT - smooth Euler Characteristic Transform
# returns approximately num_directions
compute_sect <- function(simplicial_complex, num_directions, curve_length){
  # generate equidistributed points around the sphere
  directions <- generate_equidistributed_points(num_directions)
  sect <- rep(list(NA),dim(directions)[1])
  for (i in 1:(dim(directions)[1])){
    sect[[i]] <- compute_ec_curve_direction(simplicial_complex, directions[i,], curve_length)
  }
  sect
} 

# written for the sphere; could perhaps expand this to any dimension?
# reference: https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
#How to generate equidistributed points on the surface of a sphere Markus Deserno
# use scatterplot3d to check
generate_equidistributed_points <- function(N){
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
  points
}
generate_equidistributed_points_hemisphere <- function(N){
  a <- 2*pi/N
  d <- sqrt(a)
  points <- matrix(NA,0,3) 
  
  M_theta <- round(pi/d)
  d_theta <- pi/M_theta
  d_phi <- a/d_theta
  for (i in 0:(floor(M_theta/2)-1)){
    theta <- pi*(i + 0.5)/M_theta
    M_phi <- round(2*pi*sin(theta)/d_phi)
    for (j in 0:(M_phi - 1)){
      phi <- 2*pi*j/M_phi
      point <- c( sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta))
      points <- rbind(points,point)
    }
  }
  points
}

# based on code on the EC3D.R file.
# computes the smooth euler characteristic curve in a given direction of the shape. Given a direction, in the form of a 
# unit vector, compute the vertex function by dot product; take the discrete ec curve of this simplicial complex,
# then output the smooth curve.
compute_ec_curve_direction <- function(simplicial_complex, direction, curve_length){
  vertex_function <- simplicial_complex$Vertices%*%direction
  ec_curve <- compute_discrete_ec_curve(simplicial_complex, vertex_function, curve_length)
  ec_curve <- integrate_ec_curve(ec_curve)
}

#####
# Computes the Euler Chacteristic Curve, given the simplicial complex, the values of the vertices that give rise to a 
# Sublevel set filtration, and the length of the resulting curve.
#
# This function does not integrate the EC curve -  no centering is done.
#
# Use this function to find the critical points for reconstruction.
#
# Parameter given where if the vertex indices are in the first column or given in the row indices themselves.
#####
compute_discrete_ec_curve <- function(complex, vertex_function, curve_length, first_column_index = FALSE){
    
    V = complex$Vertices
    E = complex$Edges
    F = complex$Faces
    
    if ( dim(V)[1] != dim(vertex_function)[1] ) {
      print('The size of function should be same as the number of vertices')
      return(0)
    }
    
    # In the case where the vertex_indices are given in the first column. We pick out the birth times of each vertex in every
    # edge and face of the complex.
    if (first_column_index == TRUE){
      
      edge_birth_times = sapply(1:dim(E)[1],function(ind) max(vertex_function[ which(vertex_function[,1]%in%E[ind,]), 2 ] ))
      face_birth_times = sapply(1:dim(F)[1],function(ind) max(vertex_function[ which(vertex_function[,1]%in%F[ind,]), 2 ] ))
    }else{
      
      #Apply a function on indices on this vector? Count the number of edges that are present given
      # the value of the filtration (get max vertex value at each edge, face).
      edge_birth_times = sapply(1:dim(E)[1],function(ind) max(vertex_function[ E[ind,] ]))
      face_birth_times = sapply(1:dim(F)[1],function(ind) max(vertex_function[ F[ind,] ]))
    }
    
    if (first_column_index == FALSE){
      threshold = seq(from=min(vertex_function),to=max(vertex_function),length =curve_length+1)
    } else{
      threshold = seq(from=min(vertex_function[,2]),to=max(vertex_function[,2]),length =curve_length+1)
    }
    
    
    ec_curve = matrix(0,curve_length+1,2)
    ec_curve[,1] <- threshold
    
    # count how many of each object is born given the threshold time.
    for ( i in 1:length(threshold) ) {
      v = length(which(vertex_function <= threshold[i]))
      e = length(which(edge_birth_times <= threshold[i]))
      f = length(which(face_birth_times <= threshold[i]))      
      ec_curve[i,2] <- v - e + f;
    }
    return(ec_curve)
}

#To integrate this, subtract the mean,and compute the cumulative sum
integrate_ec_curve <- function(ec_curve){
  length <- length(ec_curve[,2])
  ec_curve[,2] <- ec_curve[,2] - mean(ec_curve[,2])
  ec_curve[,2] <- cumsum(ec_curve[,2])*(ec_curve[length,1] - ec_curve[1,1])/length
  
  ec_curve
}

#### Miscellaenous Functions ####


#Computing Rotation Matrix (Want to get from your direction (x) to your desired direction (y, in our case (1,0,0))
rotation_matrix=function(x,y){
    #if (all.equal(x[2:3],c(0,0))){
    #    scalar=1/x[1]
    #    return(diag(scalar,3))
    #}
    u=x/sqrt(sum(x^2))
    v=y-sum(x*y)*u
    v=v/sqrt(sum(v^2))
    cost=sum(x*y)/sqrt(sum(x^2))/sqrt(sum(y^2))
    sint=sqrt(1-cost^2);
    rot_mat=diag(length(x)) - u %*% t(u) - v %*% t(v) + 
        cbind(u,v) %*% matrix(c(cost,-sint,sint,cost), 2) %*% t(cbind(u,v))
    return(rot_mat)
}

#### Processing OFF files from a Directory
create_ec_matrix=function(directory,direction,len){
  curve_length = len
  file_names=list.files(path=directory,full.names = TRUE)
  number_files=length(file_names)
  data <- matrix(NA,nrow=number_files,ncol = curve_length)
  for (i in 1:number_files){
    off=process_off_file_v3(file_names[i])
    vertex_function=off$Vertices%*%direction
    curve <- compute_discrete_ec_curve(off, vertex_function, len-1, first_column_index = FALSE)
    curve <- integrate_ec_curve(curve)
    data[i,]=curve[,2]
  }
  return(data)
}
create_comparison_matrix=function(directory1,directory2,direction,len){
  matrix_1=create_ec_matrix(directory1,direction,len)
  matrix_2=create_ec_matrix(directory2,direction,len)
  matrix_1=cbind(rep(1,dim(matrix_1)[1]),matrix_1)
  matrix_2=cbind(rep(0,dim(matrix_2)[1]),matrix_2)
  final_matrix=rbind(matrix_1,matrix_2)
  return(final_matrix)
}

# extend these functions to multiple directions:
create_ec_matrix_mult_d=function(directory,directions,len){
  curve_length = len
  file_names=list.files(path=directory,full.names = TRUE)
  number_files=length(file_names)
  data <- matrix(NA,nrow=number_files,ncol = curve_length*dim(directions)[1])
  for (i in 1:number_files){
    off=process_off_file_v3(file_names[i])
    curve_mult_d <- matrix(NA,nrow = 1,ncol=0)
    for (j in 1:dim(directions)[1]){
      vertex_function=off$Vertices%*%directions[j,]
      curve <- compute_discrete_ec_curve(off, vertex_function, len-1, first_column_index = FALSE)
      curve <- integrate_ec_curve(curve)
      curve_mult_d <- cbind(curve_mult_d,t(curve[,2]))
    }
    data[i,] <- curve_mult_d
  }
  return(data)
}
create_comparison_matrix_mult_d=function(directory1,directory2,directions,len){
  matrix_1=create_ec_matrix_mult_d(directory1,directions,len)
  matrix_2=create_ec_matrix_mult_d(directory2,directions,len)
  matrix_1=cbind(rep(1,dim(matrix_1)[1]),matrix_1)
  matrix_2=cbind(rep(0,dim(matrix_2)[1]),matrix_2)
  final_matrix=rbind(matrix_1,matrix_2)
  return(final_matrix)
}