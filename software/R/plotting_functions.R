#library(RColorBrewer)
#library(Rvcg)
#library(pracma)
#library(matlib)
#library(rgl)

#' Convert Mesh3d File to List
#' @description  converts a mesh3d file to a list format that's more friendly to our ECT functions.
#' @export
#' @param mesh (mesh3d): The input mesh to transform to a list.
#'
#' @return complex (list): The mesh converted to list form.
convert_off_file = function(mesh){
  vertices=as.matrix(t(mesh$vb)[,1:3])
  faces=as.matrix(t(mesh$it))
  edges=vcgGetEdge(mesh)
  edges=as.matrix(edges[,1:2])
  complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
  return(complex)
}


##########################################################################################
##########################################################################################
##########################################################################################
### Heatmap Code ###

#' Computes the birth times of the vertices on the shape.
#' @export
#' @description Given a specified number of discrete cuts and a vector of importance values, the function computes the
#' birth times of each vertex on the shape by iterating through the importance values in descending order and running the reconstruction procedure.
#' The vertex birth time is the first importance value at which the vertex is reconstructed.
#'
#'@param dir (nx3 matrix):  The matrix of directions for which the (S/D) EC curve were computed over.
#'@param complex (list) : The list containing metadata of the Vertices, Edges, and Faces of the mesh (use process_off_file_v3 to obtain).
#'@param rate_vals (vector) : Vector of variable importances for each sub-level set across each direction in a given cone.
#'@param len (int) : The number of sub-level sets which the (S/D) EC curve were computed in each direction.
#'@param cuts (int) How many steps to reconstruct when descending the importance values. A higher cut = more fine granularity on the birth times.
#'@param cone_size (int) : The number of directions in each cone.
#'@param ball_radius (float) : The radius of the ball to compute the (S/D) EC on; if you want the curve to be computed relative to the shape, don't touch this parameter.
#'@param ball (boolean) : Determining whether or not to compute the (S/D) EC curve on a ball for uniform comparisons.
#'@param radius (int) : The number of sub-level sets "before" and "after" the selected sub-level sets we want to include.
#'
#'@return vert_matrix (matrix) A matrix detailing at which importance value the vertex was reconstructed, and which "cut" the vertex was reconstructed relative to the other vertices.
reconstruct_vertices_on_shape = function(dir, complex, rate_vals, len, cuts=10, cone_size, ball_radius,ball = TRUE, radius = 0){
  vert_matrix = matrix(0,nrow = dim(complex$Vertices)[1], ncol = 2)
  cut = cuts
  reconstructed_vertices = c()
  for (threshold in quantile(rate_vals,probs = seq(1,0,length.out = cuts)) ){
      selected_vertices = compute_selected_vertices_cones(dir = dir, complex = complex, rate_vals = rate_vals, len = len, threshold = threshold,
                                                          cone_size = cone_size,ball_radius = ball_radius, ball = ball, radius = radius)
      selected_vertices = setdiff(selected_vertices,reconstructed_vertices)
      vert_matrix[selected_vertices,1] = cut
      vert_matrix[selected_vertices,2] = threshold
      cut = cut-1
      reconstructed_vertices = c(reconstructed_vertices,selected_vertices)
      if (length(reconstructed_vertices) == dim(complex$Vertices)[1]){
        break
      }
}
  return(vert_matrix)
}

#' Computes the  birth times of the faces on the shape.
#' @export
#' @description Given a specified number of discrete cuts and a vector of importance values, the function computes the
#' birth times of each face on the shape by iterating through the importance values in descending order and running the reconstruction procedure.
#' The face birth time is the first importance value at which a vertex corresponding to a face face is reconstructed.
#'
#'@param dir (nx3 matrix):  The matrix of directions for which the (S/D) EC curve were computed over.
#'@param complex (list) : The list containing metadata of the Vertices, Edges, and Faces of the mesh (use process_off_file_v3 to obtain).
#'@param rate_vals (vector) : Vector of variable importances for each sub-level set across each direction in a given cone.
#'@param len (int) : The number of sub-level sets which the (S/D) EC curve were computed in each direction.
#'@param cuts (int) : How many steps to reconstruct when descending the importance values. A higher cut = more fine granularity on the birth times.
#'@param cone_size (int) : The number of directions in each cone.
#'@param ball_radius (float) : The radius of the ball to compute the (S/D) EC on; if you want the curve to be computed relative to the shape, don't touch this parameter.
#'@param ball (boolean) : Determining whether or not to compute the (S/D) EC curve on a ball for uniform comparisons.
#'@param radius (int) : The number of sub-level sets "before" and "after" the selected sub-level sets we want to include.
#'
#'@return face_matrix (matrix) A matrix detailing at which importance value the face was reconstructed, and which "cut" the face was reconstructed relative to the other faces.
reconstruct_faces_on_shape = function(dir, complex, rate_vals, len, cuts=10, cone_size, ball_radius,ball = TRUE, radius = 0){
  face_matrix = matrix(0,nrow = dim(complex$Faces)[1], ncol = 2)
  cut = cuts
  reconstructed_faces = c()
  for (threshold in quantile(rate_vals,probs = seq(1,0,length.out = cuts)) ){
      selected_faces = compute_selected_faces_cones(dir = dir, complex = complex, rate_vals = rate_vals, len = len, threshold = threshold,
                                                          cone_size = cone_size,ball_radius = ball_radius, ball = ball, radius = radius)
      selected_faces = setdiff(selected_faces,reconstructed_faces)
      face_matrix[selected_faces,1] = cut
      face_matrix[selected_faces,2] = threshold
      cut = cut-1
      reconstructed_faces = c(reconstructed_faces,selected_faces)
      if (length(reconstructed_faces) == dim(complex$Faces)[1]){
        break
      }
}
  return(face_matrix)
}

#' Reconstructs and finds vertex colors for a set of meshes.
#'
#' @export
#' @description Given a set of meshes, color function, and importance values, number of cuts to reconstruct on,the function finds the birth times of all the
#' vertices and finds vertex colors using the color function provided.
#'
#'@param dir_name (string): The input directory containing the meshes to loop over
#'@param cuts (int): How many steps to reconstruct when descending the importance values. A higher cut = more fine granularity on the birth times.
#'@param pset (list): A list containing the metadata about the analyses (directions, number of sub level sets, directions per cone)
#'@param comp (list): A list of elements containing metadata about the variable selection (primarily variable importances)
#'@param colfunc (function): The function for the color gradients.
#'
#'@return mesh_list (list) A list containing the colors of the meshes derived from the iterative reconstruction procedure.
get_heat_colors = function(dir_name, cuts, pset, comp, colfunc){
  mesh_list = list()
  for (k in 1:length(list.files(dir_name,full.names = TRUE))){
    file_name = list.files(dir_name,full.names = TRUE)[k]
    print(paste('On File', file_name))
    mesh = process_off_file_v3(file_name)
    heat = reconstruct_vertices_on_shape(dir = pset$dirs,complex = mesh,rate_vals = comp$Rate2[,2],
                                                len = pset$len,cuts = cuts,cone_size = pset$directions_per_cone,ball_radius = ball_radius, ball = ball,radius =1)
    heat_colors = colfunc(1 + max(heat[,1]) - min(heat[,1]))[1 + heat[,1] - min(heat[,1])]
    mesh_list[[k]] = heat_colors
  }
  return(mesh_list)

}

#' Creates a color bar:
#'
#' @export
#' @description Taken from http://www.colbyimaging.com/wiki/statistics/color-bars. Creates a color bar from a given color scale and range of tick values.
#'
#' @param lut (function) : The color scale.
#' @param min (float) : Minimum value on the color scale
#' @param max (float): Maximum value of the color scale
#' @param nticks (int) : Number of ticks to have on the color bar.
#' @param ticks (vector): The vector of ticks (the user can pass in a vector of custom axes ticks)
#' @param title (string): The title of the color bar.
#'
#' @return color_bar (figure) : A color bar.
color.bar <- function(lut, min, max=-min, nticks=11, ticks=seq(min, max, len=nticks), title='') {
  scale = (length(lut)-1)/(max-min)

  dev.new(width=1, height=5)
  par(mar=c(5,6,4,1)+.1)
  plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
  axis(2, ticks, las=1,cex.axis = 1.4)
  for (i in 1:(length(lut)-1)) {
    y = (i-1)/scale + min
    rect(0,y,10,y+1/scale, col=lut[i], border=NA)
  }
}
#------------------------------------------------------------------------------------------------------------------------------#
# map rate_values to vertices; output is a list.
vertex_to_feature_dict <- function(complex, dirs, curve_length, ball_radius){
  num_vertices <- dim(complex$Vertices)[1]
  f <- functional::Curry(vertex_to_projections, dirs = dirs, curve_length = curve_length, ball_radius = ball_radius)

  t(apply(complex$Vertices, 1, f))
}

# maps a vertex to the projections in all directions
vertex_to_projections <- function(vertex, dirs, curve_length,ball_radius){
  f <- functional::Curry(vertex_to_dir_projection, vertex = vertex, dirs = dirs,
                         curve_length = curve_length, ball_radius = ball_radius)

  purrr::map_int(1:(dim(dirs)[1]), f)
}

# idx represents the index in the set of directions
# returns for each vertex in the complex, the index of the projection in the given direction
vertex_to_dir_projection <- function(vertex, idx, dirs, curve_length, ball_radius){
  vtx_projection <- vertex%*%dirs[idx,]
  buckets <- seq(-ball_radius,ball_radius,length.out = curve_length+1)

  # map vertex projection to the feature index
  projection_bucket <- cut(vtx_projection, buckets, labels = FALSE)

  # update index to reflect rate values
  projection_bucket <- as.integer(projection_bucket + (idx - 1)*curve_length)

  projection_bucket
}


# color the vertices of the complex based on how many RATE-selected projections it lies in
get_projection_heatmap <- function(complex, directions, rate_values,
                                   curve_length,ball_radius, threshold = 1/length(rate_values)){
  rate_selected_features <- which(rate_values > threshold)
  vertex_to_feature <- vertex_to_feature_dict(complex, directions, curve_length, ball_radius)

  # each entry in the list is
  vertex_to_num_projections <- matrix(0,ncol = dim(complex$Vertices)[1])
  for (i in 1:dim(complex$Vertices)[1]){
    vertex_to_num_projections[i] = sum(vertex_to_feature[i,] %in% rate_selected_features)
  }
  vertex_to_num_projections
}


# color the vertices of the complex based on sum of RATE values on each projection it lies in
get_RATE_weighted_heatmap <- function(complex, directions, rate_values,
                                      curve_length,ball_radius, threshold = 1/length(rate_values)){
  rate_selected_features <- which(rate_values > threshold)
  vertex_to_feature <- vertex_to_feature_dict(complex, directions, curve_length, ball_radius)

  # each entry in the list is
  vertex_to_rate_weighted_projections <- matrix(0,ncol = dim(complex$Vertices)[1])
  for (i in 1:dim(complex$Vertices)[1]){
    vertex_to_rate_weighted_projections[i] = sum(rate_values[which(vertex_to_feature[i,] %in% rate_selected_features)])
  }
  vertex_to_rate_weighted_projections
}

# Wrapper for plotting teeth within each group #

plot_results_teeth_simple=function(files, features1, features2, color1, color2, alpha1=0.65, alpha2=0.65,
                                   dir1, dir2, len=65, level=20, slices=25, n=10, thresh = 1.00,
                                   directions_per_cone = 5){
  #Loop through files in the directory
  for (i in 1:length(files)){
    #We set direction as (0,0.75,1)
    file1=vcgImport(files[i])
    file_1=process_off_file_v3(files[i])
    #Try to see if this file has any critical points, if it doesn't just plot the tooth.
    vert1 = compute_selected_vertices_cones(dir = dir1, complex =file_1, rate_vals = features1, len = len, threshold = (thresh/(length(features1))),
                                            cone_size = directions_per_cone,ball_radius = ball_radius, ball = FALSE)
    vert2 =  compute_selected_vertices_cones(dir = dir2, complex =file_1, rate_vals = features2, len = len, threshold = (thresh/(length(features2))),
                                             cone_size = directions_per_cone,ball_radius = ball_radius, ball = FALSE)
    intersected = intersect(vert1,vert2)
    fc3 <- colorRampPalette(c(color1,color2))
    colors = rep('white', dim(veg1$vb)[2])
    colors[setdiff(vert1,vert2)] =fc3(10)[2]
    colors[setdiff(vert2,vert1)] =fc3(10)[9]
    colors[intersected] = fc3(10)[6]
    plot3d(file1, colors = colors,axes = FALSE, xlab = '',ylab = '', zlab = '')
    #Rotate the tooth for a view of the 'bottom'
    rotation_matrix=matrix(c(0.99972576,0.02127766,0.00978078,0,0.01589701,-0.92330807,0.38373107,0,0.01719557,-0.38347024,-0.92339307,0,0,0,0,1),ncol=4,byrow=TRUE)
    rgl.viewpoint(userMatrix = rotation_matrix)
  }
}

