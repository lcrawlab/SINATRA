library(Rvcg)
library(rgl)
library(FNN)
library(pracma)
library(matlib)
library(far)


#' Generates 4 directions in a small cone around a specified direction.
#'
#' \code{generate_perturbed_directions} returns a list of the input direction and the four directions.
#'
#' The function first computes two orthogonal vectors to the input direction, then perturbs the original direction
#' by rotating the original direction an angle of phi and -phi in the direction of the two orthogonal vectors, which
#' generates four additional directions in a cone around the original.
#'
#'
#' @param direction A vector. Should ideally be normalized so norm is one.
#' @param phi A scalar. The angle parameter for the angle between the given direction and its adjacent directions.
#' @return The output is a list of length 5, with direction is one entry of the list.
#'
generate_perturbed_directions=function(direction,phi=0.001){
  # initialize the list
  directions=list()

  # generate orthonormal directions, via Gram Schmidt
  basis=orthonormalization(direction,norm=TRUE)
  v0=basis[,1]
  v1=basis[,2]
  v2=basis[,3]
  isna=unique(v2)
  if (length(isna)==1){
    basis=orthonormalization(cbind(v0,v1))
    v0=basis[,1]
    v1=basis[,2]
    v2=basis[,3]
  }
  c1=cos(phi)
  c2=sin(phi)

  # Generate additional vectors in the direction of v1
  direction2=c1*v0+c2*v1
  direction3=c1*v0-c2*v1

  # Generate additional vectors in the direction of v2
  direction4=c1*v0+c2*v2
  direction5=c1*v0-c2*v2

  # generate list
  directions[[1]]=direction
  directions[[2]]=direction2
  directions[[3]]=direction3
  directions[[4]]=direction4
  directions[[5]]=direction5
  return(directions)
}




#' Returns the critical points projections given a direction and the simplicial complex.
#'
#' \code{find_critical_points_curve} Returns a list of the values at which the ec curve changes, for the directions in a cone around the input direction.
#'
#' The function is to compute the critical point projections for a cone of directions.
#' Can be used more generally so long as the input is a list. For each direction, the function computes the unsmoothed Euler Characteristic
#' Curve and computes the indices and corresponding coordinate for which the EC changes.
#'
#' @param directions A list. Each entry is a separate direction (a unit vector), in matrix form. This is the output of
#'  \code{generate_perturbed_directions}.
#' @param off A simplicial complex. Can be obtained by calling methods from SECT on an OFF file.
#' @param len Integer. Length of the ec curve.
#' @return The output is a list: each entry i consists of the projections of critical points onto direction i, defined as
#' the index where the ec curve in direction i changes value. The coordinate of the projection (the distance from the origin
#' to the projected) is also recorded.
find_critical_points_curve=function(direction, off, steps,phi){
  critical_points=list()

  # Generate a cone of directions around
  directions=generate_perturbed_directions(direction = direction,phi=phi)

  # For each
  for (i in 1:5){
    direction=directions[[i]]
    vertex_function=off$Vertices%*%direction
    ec_curve=compute_discrete_ec_curve(off,vertex_function,steps)
    num_t=dim(ec_curve)[1]
    crit_points=c()
    interval=c()

    # find the indices at which the EC changes
    for (j in 1:(num_t-1)){
      current_value=ec_curve[j,2]
      next_value=ec_curve[(j+1),2]
      if (next_value!=current_value){
        crit_points=c(crit_points,ec_curve[(j+1),1])
        interval=c(interval,j+1)}

    }

    # Get the critical points
    if (length(crit_points)==0){
      #print('No Critical Points!')
    }
    else{
      crit_points=as.matrix(cbind(interval,crit_points))
      #print(crit_points)
      critical_points[[i]]=crit_points
    }
  }
  return(critical_points)
}






#' Critical Point Finding Algorithm. Returns the coordinates of critical points given a direction and the simplicial complex.
#'
#' \code{find_critical_point_coordinates} computes the coordinates of critical points given a direction and the desired simplicial
#' complex. The number of steps denotes the number of sublevel sets to consider when computing the EC curve, and phi is the initial
#' choice of angle for the cone of directions around the input.
#'
#' This function computes the coordinates of critical points of the desired simplicial complex in a specified direction; that is, the
#' points for which the EC changes in this direction. There is also some further interpretation in Morse Theory to be further explained
#' in the paper. Given a direction, we generate a cone of 5 directions around and including the input direction. For each of these
#' directions, we compute the critical point projections and the indices of the ec curves at which they occur. We then shrink phi until
#' each of these directions outputs the same number of critical points with the assumption that these projections refer to the same points.
#' Once this happens, we solve the overdetermined putative critical point equation to back out the coordinates of the critical points.
#' See the paper for a discussion about the assumptions.
#'
#' KNOWN ISSUES: this function returns critical points that lie on a straight line. We think it may be due to sensitivity to the parameters,
#' especially the number of sublevel sets taken for the EC curve, or perhaps numerical sensitivity when the value of phi is extremely small,
#' which occurs when this function is run. The issue could also be with solving the putative critical point equation, which could algebraically
#' explain why the critical points lie on a subspace. It is currently unclear how the straight line is related to the matrix in the putative
#' critical point equation. The current fix is to take out the returned critical points that lie too far away from the shape. Another idea for
#' why this issue occurs is perhaps we are working with the wrong coordinates. In particular the values in each direction at which the EC curve
#' changes is not in cartesian standard coordinates.
#'
#' @param direction a vector, in matrix type.
#' @param off a simplicial complex: a list containing vertices, faces, and edges. Obtained by calling methods from the SECT_Functions file.
#' @param steps Integer. Length of the ec curve.
#' @param phi Float. Degree of the angle between input direction and the original.
#' @param rejection_sampling(TRUE/FALSE) A boolean that determines if the rejection sampling procedure will happen. (Take out the critical points far away from the shape).
#' @return A matrix. The matrix is n x 4 where n is the number of critical points found. The last three columns are the x,y,z coordinates of
#' the critical points. The first column is the index in the ec curve that corresponds to that critical point.
find_critical_point_coordinates=function(direction,off,steps,phi,rejection_sampling=TRUE){

  # what is this block used for?
  directions=generate_perturbed_directions(direction = direction,phi=phi)
  basis=orthonormalization(direction,norm=TRUE)
  v0=basis[,1]
  v1=basis[,2]
  v2=basis[,3]
  isna=unique(v2)
  if (length(isna)==1){
    basis=orthonormalization(cbind(v0,v1))
    v0=basis[,1]
    v1=basis[,2]
    v2=basis[,3]
  }
  #inverse_basis=solve(basis)

  #Find all the critical points projections for each direction in a cone around the input (where the curve changes direction)
  crit_points=try(find_critical_points_curve(direction=direction,off=off,steps=steps,phi=phi))

  # Next, assert that each set of critical points per direction is equal (to ensure phi is small enough). Recall that this is for the five
  # directions in the cone around the input.
  if (length(crit_points)==0){
    #print('No Critical Points!')
    return(0)
  }
  if (length(crit_points)!=5){
    phi=phi*0.1
    #print('Phi not small enough, multiplying by 0.1')
    return(find_critical_point_coordinates(direction=direction,off=off,steps=steps,phi = phi))
  }
  benchmark_length=dim(crit_points[[1]])[1]
  critical_point_lengths=c()
  direction_matrix=matrix(NA,ncol=3,nrow=5)
  for (i in 1:5){
    critical_point_lengths=c(critical_point_lengths,dim(crit_points[[i]])[1])
    direction_matrix[i,]=directions[[i]]
  }
  test_equal_length=(benchmark_length==critical_point_lengths)
  equal_lengths=unique(test_equal_length)


  #Once the same number of critical points are captured in each direction, solve the putative critical point equation; otherwise shrink phi
  if (length(equal_lengths)==1){
    if (equal_lengths==TRUE){
      intervals=crit_points[[1]][,1]
      xs=c()
      ys=c()
      zs=c()
      lens=dim(crit_points[[1]])[1]

      # For each of the critical points, compute their coordinates.
      for (i in 1:lens){
        paired_crit_points=c()
        for (k in 1:5){
          paired_crit_points=c(paired_crit_points,crit_points[[k]][i,2])
        }

        # solve the critical point equation
        solution=try(gaussianElimination(A = direction_matrix, B=as.matrix(paired_crit_points)))
        if(inherits(solution,'try-error')){
          #print(direction_matrix)
          #print(matrix(paired_crit_points))
          solution=gaussianElimination(A = direction_matrix[1:length(paired_crit_points),], B=as.matrix(paired_crit_points))
          xs=c(xs,solution[1,4])
          ys=c(ys,solution[2,4])
          zs=c(zs,solution[3,4])}
        else{
          xs=c(xs,solution[1,4])
          ys=c(ys,solution[2,4])
          zs=c(zs,solution[3,4])
        }

      }
      total_crit_coord=cbind(intervals,xs,ys,zs)
      remove=c()

      # Current fix to the issue of critical points lying on a line: take out the results not near the shape.
      if (rejection_sampling==TRUE){
        for (j in 1:dim(total_crit_coord)[1]){
          dist=apply(X = off$Vertices,MARGIN = 1,FUN = difference,y=total_crit_coord[j,2:4])
          if (min(dist)>0.1){
            remove=c(remove,j)
            #print(min(dist))
          }
        }
      }
      if (length(remove)!=0){
        total_crit_coord=total_crit_coord[-remove,]
      }
      return(total_crit_coord)
    }
    else{
      phi=phi*0.1
      #print('Phi not small enough, multiplying by 0.1')
      return(find_critical_point_coordinates(direction=direction,off=off,steps=steps,phi = phi))
    }
  }
  else{
    phi=phi*0.1
    #print('Phi not small enough, multiplying by 0.1')
    return(find_critical_point_coordinates(direction=direction,off=off,steps=steps,phi = phi))
  }
}





#' Finds critical points for multiple directions.
#'
#' \code{find_critical_points_multiple_directions} computes the coordinates of critical points given a list of directions and the desired simplicial
#' complex. The number of steps denotes the number of sublevel sets to consider when computing the EC curve, and phi is the initial angle.
#'
#'
#' @param dir A matrix of directions, n x 3. Each row contains a direction for which to compute the critical points.
#' @param off a simplicial complex: a list containing vertices, faces, and edges. Obtained by calling methods from the SECT_Functions file.
#' @param steps Integer. Length of the ec curve.
#' @param phi Float. Degree of the angle between input direction and the original.
#' @param rejetion_sampling (True/FALSE) A boolean that determines if we will do rejection sampling or not.
#' @return A list. Each entry corresponds to a different direction, and the entry is an n x 4 matrix, where n is the number of critical points found.
#'  The last three columns are the x,y,z coordinates of the critical points. The first column is the index in the ec curve that corresponds to that
#'  critical point; this is adjusted so that it matches the direction that the critical point came from, giving it a unique identifier in terms of which
#'  direction and where on the ec curve that critical point corresponds to.
find_critical_points_multiple_directions=function(dir,off,steps,phi,rejection_sampling=TRUE){
  critical_points=list()
  total_directions=dim(dir)[1]
  for (i in 1:total_directions){
    direction=dir[i,]
    crit_points=find_critical_point_coordinates(direction = direction,off = off,steps = steps,phi = phi,rejection_sampling = rejection_sampling)
    if (length(crit_points)>2){
      crit_points=matrix(crit_points,ncol=4)
      crit_points[,1]=crit_points[,1]+(steps*(i-1))
    }
    critical_points[[i]]=crit_points
  }
  return(critical_points)
}




#' Extracts coordinates of critical points from a list.
#'
#' \code{get_all_critical_points} is a simple wrapper for \code{find_critical_points_multiple_directions}; it takes the output of the latter function
#' and extracts the coordinates of the critical points for each direction.
#'
#' @param critical_point_list A list. Each entry is a matrix of critical points, the first column are the indices of the ec curve and the remaining
#' columns are the coordinates of the critical points.
#'
#' @return A matrix of solely the coordinates of the critical points.
get_all_critical_points=function(critical_point_list){
  critical_points=c(0,0,0,0)
  for (i in 1:length(critical_point_list)){
    critical_points=rbind(critical_points,critical_point_list[[i]])
  }

  if (length(critical_point_list) > 1){
    critical_points <- unique(critical_points[-1,-1])
  } else{
    critical_points <- critical_points[-1,-1]
  }
  critical_points=matrix(critical_points,ncol=3)
  return(critical_points)
}





#' Extracts coordinates of critical points from a list.
#'
#' \code{get_all_critical_points} is a simple wrapper for \code{find_critical_points_multiple_directions}; it takes the output of the latter function
#' and extracts the coordinates of the critical points for each direction along with the indices.
#'
#' @param critical_point_list A list. Each entry is a matrix of critical points, the first column are the indices of the ec curve and the remaining
#' columns are the coordinates of the critical points.
#'
#' @return A matrix of the coordinates of the critical points along with the indices.
get_all_critical_points_feature_selection=function(critical_point_list){
  critical_points=c(0,0,0,0)
  for (i in 1:length(critical_point_list)){
    critical_points=rbind(critical_points,critical_point_list[[i]])
  }

  # Take the unique critical points
  if (length(critical_point_list) > 1){
    critical_points <- unique(critical_points[-1,])
  } else{
    critical_points <- critical_points[-1,]
  }

  # what does this do?
  critical_points=matrix(critical_points,ncol=4)
  critical_points=critical_points[critical_points[,1]!=0,]
  return(critical_points)
}




#' The Reconstruction Procedure, given the indices deemed important.
#'
#' \code{critical_points_feature_selection} returns the coordinates of the critical points that are chosen by the feature selection procedure.
#'
#' The code first checks if the number of directions in the critical points input and directions input match up. If so, we pick out all the
#' critical points that have indices contained within the desired indices, the output of the feature selection algorithm. Essentially, each critical
#' point comes with an identifier unique to the direction and index on that direction's ec curve. If that identifier matches with the index chosen
#' by the feature selection algorithm, we select that critical point. We take a union of sorts: if a critical point has a projection onto a direction
#' that is chosen by feature selection, we select that critical point.
#'
#' This procedure can be changed, but at least should be discussed; we could take the intersection for example, but there is no guarantee that
#' any critical point satisfies this property. Could revive the evidence threshold idea.
#'
#' @param critical_points (A list of length n). Each entry are the critical points for each direction in matrix form, with indices in the first column.
#' @param indices (A list/matrix). The output of the feature selection algorithm. The ec curves in each direction are concatenated into a long vector, and feature
#' selection algo selects indices that are important for each direction.
#' @param direction (A n x 3 matrix). A matrix of directions; should be the ones that generated the input critical points as well.
#'
#' @return Returns the coordinates of the selected critical points.
critical_points_feature_selection=function(critical_points,indices,direction){
  num_dir=dim(direction)[1]
  num_crit_points=length(critical_points)
  crit_points_total=c(0,0,0,0)

  # check for the correct number of directions
  if (num_dir!=num_crit_points){
    print('Error: Collection of Critical Points Does not Match with Number of Directions')
    print(num_dir)
    print(num_crit_points)
    return(0)
  }
  # otherwise, choose the critical points by the union procedure described in the function doc.
  else{
    crit_points=get_all_critical_points_feature_selection(critical_points)
    crit_points_total=crit_points[crit_points[,1] %in% indices, ]
  }

  # out of the chosen critical points, take the unique ones.
  if (is.null(dim(crit_points_total)[1])){
    return(vector(mode="numeric", length=0))
  }
  else{
    crit_points_total=crit_points_total[,-1]
    crit_points_total=unique(crit_points_total)
  }
  return(crit_points_total)
}





#### Functions for Plotting ####
#' The plotting functions: (create a convex hull of the vertices closest to the critical points)
#'
#' \code{create_hull} Plot the Convex Hull of the points given some vertices of the convex hull.
#'
#' The code transverses through the matrix of critical points, finds the 'n' closest vertices (n is a parameter) for each critical point.
#' For each point, the code grabs the coordinates of the 'n' closest vertices, and perturbs them slightly.
#' The perturbation is done becase the rgl plotting procedure doesn't allow for overlapping shapes on the exact same coordinates.
#' The perturbed coordinates of the vertices are then created into a convex hull, and then plotted.
#' If necessary, the critical points can be plotted too.
#'

#'
#' @param points A (nx3 matrix). Each entry are the critical points in matrix form. These points are what we create the convex hull around.
#' @param vertices A (mx3 matrix). The vertices of the simplicial complex for which the function creates a convex hull around.
#' @param n (integer) A parameter to choose, how big the convex hull is (how many of the n closest vertices).
#' @param color (string) Color of the convex hull.
#' @param alpha (float (0-1)): The intensity of the color of the convex hull
#' @param plot_points (TRUE/FALSE): A boolean to see if the user wants to plot the critical points themselves, in addition to the hull.
#'
#' @return Plots the convex hull of the critical points in the RGL interface.
create_hull=function(points,vertices,n,color='blue',alpha=0.3,plot_points=FALSE){
  number_crit_points=dim(points)[1]
  for (i in 1:number_crit_points){
    dummy_point=points[i,]
    dummy_point=matrix(dummy_point,ncol=3)
    closest_points_index=knnx.index(vertices,dummy_point,k=n)
    closest_points=vertices[closest_points_index,]
    perturb=matrix(rnorm(n*3,mean=0,sd=0.001), n, 3)
    closest_points=closest_points+perturb
    ts.surf=t(convhulln(closest_points))
    #To make sure overlapping areas get covered
    rgl.triangles(closest_points[ts.surf,1],closest_points[ts.surf,2],closest_points[ts.surf,3],col=color,alpha=alpha)
    rgl.points(closest_points[1,],color=color,alpha=1,size=5)
    if (plot_points==TRUE){
      rgl.points(dummy_point,color=color,size=10)
    }
  }
}

#' \code{plot_points} Plot the critical points .
#'
#' The code plots the critical points.
#'
#'
#' @param points A (nx3 matrix). Each entry are the critical points in matrix form. These points are what we create the convex hull around.
#' @param point_size (int): Specifies the point size of the critical points.
#' @param point_color (string) Color of the shapes.
#' @param alpha (float (0-1)): The intensity of the color of the points.
#'
#' @return Plots the Critical points.
plot_points=function(points,point_size=5,point_color='red'){
  rgl.points(points)
  rgl.points(points,col=point_color,size=point_size)
}

#' \code{plot_complex_with_points} Plot the 3d shape, and the critical points (convex hull or points themselves).
#'
#' The code plots the 3d shape, and the plotting of the critical points. There are a few plotting parameters, determining the shading and transparency of the shapes.
#' The critical points provided are then plotted too. They can be plotted either as points, or as a convex hull. These are given as boolean conditions.
#'

#'
#' @param points A (nx3 matrix). Each entry are the critical points in matrix form. These points are what we create the convex hull around.
#' @param n (integer) A parameter to choose, how big the convex hull is (how many of the n closest vertices), if chosen.
#' @param main_color (string) Color of the shape.
#' @param crit_color (string) Color of the convex hull/shape.
#' @param alpha_two (float (0-1)): The intensity of the color of the shape.
#' @param alpha_two (float (0-1)): The intensity of the color of the convex hull (if chosen)
#' @param off_file (mesh3d): A mesh3d file, read in from Rvcg. This is the shape that's plotted.
#' @param axes (TRUE/FALSE): A boolean to determine if the axes will be shown.
#' @param labels (TRUE/FALSE): A boolean to determine if the axes labels will be shown.
#' @param plot_points (TRUE/FALSE): A boolean to see if the user wants to plot the critical points themselves, in addition to the hull.
#' @param hull (TRUE/FALSE): A boolean to determine if the user wants to display the convex hull.
#' @param point_size (int): If the critical points are plotted, this specifies the point size.
#'
#' @return Plots the 3d shape,and the plotting of the critical points/or the convex hull of the critical points.
plot_complex_with_points=function(points,vertices,n=100,main_color='white',crit_color='blue',alpha_one=0.5,alpha_two=0.3,off_file,axes=FALSE,labels=FALSE,plot_points=FALSE,hull=TRUE,point_size=10){
  colors=rep(main_color,dim(off_file$vb)[2])
  vertices=t(off_file$vb)[,1:3]
  plot3d(off_file,xlab='',ylab='',zlab='',type=c('shade'),col=colors,alpha=alpha_one,axes=axes,labels=labels)
  aspect3d(1,1,1)
  if (hull==TRUE){
    if (length(points)!=0){
      create_hull(points=points,vertices,n=n,color=crit_color,alpha=alpha_two,plot_points=plot_points)
    }
  }
  else{
    plot_points(points,point_size=point_size,point_color=crit_color)
  }
}






##### Old / Misc Reconstruction Functions #####


#### Helper Functions ####
#' \code{norm_vec} Computes the norm of a vector .

#'
#' @param x (vector). A vector in \R^n.
#' @return The norm of the vector
norm_vec = function(x){
  scale=(sqrt(sum(x^2)))
  return(x/scale)
}
#' \code{difference} Computes the Euclidean distance between two vectors.

#'
#' @param x (vector). A vector in \R^n.
#' @param y (vector). Another vector in \R^n.
#' @return The Euclidean distance between the two vectors.
difference=function(x,y){
  sumxy=sum((x-y)^2)
  return(sqrt(sumxy))
}



create_hull_kde=function(points,vertices,n,color='blue',alpha=0.3,plot_points=FALSE){
  number_crit_points=dim(points)[1]
  for (i in 1:number_crit_points){
    dummy_point=points[i,]
    dummy_point=matrix(dummy_point,ncol=2)
    closest_points_index=knnx.index(vertices[,1:2],dummy_point,k=n)
    closest_points=vertices[closest_points_index,]
    ts.surf=t(convhulln(closest_points))
    rgl.triangles(closest_points[ts.surf,1],closest_points[ts.surf,2],closest_points[ts.surf,3],col=color,alpha=alpha)
    if (plot_points==TRUE){
      rgl.points(dummy_point,color=color,alpha=1,size=10)
    }
  }
}


plot_complex_with_points_simple=function(points,vertices,n=5,main_color='white',vert_colors,alpha_one=0.5,alpha_two=0.3,off_file,axes=FALSE,labels=FALSE,plot_points=FALSE,hull=TRUE,point_size=5){
  colors=rep(main_color,dim(off_file$vb)[2])
  vertices=t(off_file$vb)[,1:3]
  plot3d(off_file,xlab='',ylab='',zlab='',type=c('shade'),col=colors,alpha=alpha_one,axes=axes,labels=labels)
  aspect3d(1,1,1)
  if (hull==TRUE){
    if (length(points)!=0){
      create_hull_on_shape(points=points,vertices,n=n,colors=vert_colors,alpha = alpha_two,plot_points=plot_points)
    }
  }
  else{
    plot_points(points,point_size=point_size,point_color=crit_colors)
  }
}

find_evidence_colors=function(dir,complex,indices,len){

  evidence <- matrix(0,nrow = dim(complex$Vertices)[1], ncol = dim(dir)[1])

  # Count how many projections are selected for
  for(i in 1:dim(dir)[1]){
    projections <- complex$Vertices[,1:3]%*%dir[i,]
    buckets <- seq(min(projections),max(projections),length.out = len)

    #bucket these projections into curve_length number of groups; could have also solved this with the cut function
    step_length <- (max(projections) - min(projections))/len
    projection_buckets <- apply((projections - min(projections))/step_length,1, function(float) as.integer(float)) + len*(i-1)
    evidence[,i] <- sapply(projection_buckets %in% indices,as.numeric)
  }

  vertex_evidence <- rowSums(evidence)
  #cols=rainbow(len(unique(vertex_evidence)))[vertex_evidence]

  colfunc <- colorRampPalette(c("white", 'white','white',"blue"))
  colors=colfunc(max(vertex_evidence))[vertex_evidence]
  plot(rep(1,max(vertex_evidence)),col=colfunc(max(vertex_evidence)),pch=19,cex=3)
  return(colors)

}








