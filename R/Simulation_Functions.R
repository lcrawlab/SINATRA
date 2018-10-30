library("ks")
library("mvtnorm")
library('pdist')
library(MASS)
library(ks)
library(truncnorm)
library(spatstat)

#Kernel Function
rbf_kernel_func=function(X,Y,length_scale = 1){
    #Inputs: X (nxm matrix): First array of position coordinates
    #        Y (nxm matrix): Second Array of Position Coordinates
    #        length_scale (int): The constant we divide the distances by in the computation of the kernel.
    #Outputs: kernel: The Kernel/covariance matrix.
    
    D = as.matrix(pdist(X,Y))
    kernel = exp(-D/length_scale^2)
    return(kernel)
}
#Generating a collection of points with a small distance of training points.
#This is optional; it turns out simulations work better if we don't use this. 
generate_clustered_points_with_means=function(means,n=5,eta=0){
    data=c(0,0,0)
    n_clusters=dim(means)[1]
    for (i in 1:n_clusters){
        x=means[i,]
        data=rbind(data,x)
        for (j in 1:n){
            x_new=x+rnorm(3,eta,0.02)
            data=rbind(data,x_new)
        }
    }
    data=as.matrix(data)
    data=data[2:(n*n_clusters),]
    return(data)
}
#Generating GP with Training Points
generate_gp=function(grid_length=20,new_data,length_scale=5){
    #Inputs: grid_length (int) : How many points we want to generate a GP on within [0,1],
    #        new_data (nxm matrix): The training points to generate GP on.
    #        length_scale (int) The constant we divide the distances by in the computation of the kernel.
    #Outputs: posterior_matrix: a matrix of size grid_length*grid_length, with the entries as computed by the GP.
    x <- seq(0,1,length=grid_length)
    X <- expand.grid(x, x)
    # compute squared exponential kernel on pairwise values
    X <- as.matrix(X)
    new_Cov=rbf_kernel_func(X,X,length_scale)
    f_prior= mvrnorm(n=1,mu=rep(0,dim(new_Cov)[1]), Sigma=new_Cov, tol=0.1)
    prior_plot=data.frame(x1=X[,1],x2=X[,2],y = f_prior)
    prior.matrix= matrix(prior_plot[,3],nrow = grid_length, byrow = TRUE)
    #Use these Test Points
    obs_x=new_data[,1]
    obs_y=new_data[,2]
    obs_z=new_data[,3]
    new_data=cbind(obs_x,obs_y)
    pretty_new_data=cbind(new_data,obs_z)
    #Kernel for original data
    kernel_xsxs=rbf_kernel_func(X,X,length_scale)
    kernel_xsx=rbf_kernel_func(X,new_data,length_scale)
    kernel_xxs=rbf_kernel_func(new_data,X,length_scale)
    #Kernel for new data
    kernel_xx=rbf_kernel_func(new_data,new_data,length_scale)
    #The new "F"
    f.star.bar=kernel_xsx%*%solve(kernel_xx)%*%obs_z
    #New Covariance Matrix
    cov.f.star <- kernel_xsxs - kernel_xsx%*%solve(kernel_xx)%*%kernel_xxs
    #Generate Posterior Samples
    #f_posterior= t(rmvnorm(1,f.star.bar, cov.f.star, method = "eigen"))
    f_posterior= mvrnorm(n=1,mu=f.star.bar, Sigma=cov.f.star, tol=0.1)
    posterior=data.frame(x1=X[,1],x2=X[,2],y = f_posterior)
    #yy1.matrix= matrix(t(Y)[,1],nrow = 20, byrow = TRUE)
    posterior_matrix= matrix(posterior[,3],nrow = grid_length, byrow = TRUE)
    return(posterior_matrix)
}
#Generating multiple GP's.
create_data_gp <- function(num_sim, data_1,data_2, grid_size=20,length_scale=5,dir,curve_length = 50){
    #Inputs:  num_sim: observations for each class
    #         data_1 (nxm array): the 'causal' points for class 1
    #         data_2 (nxm array) the 'causal' points for class 2
    #         grid_size (int) : number of training points on [0,1]
    #         length_scale (int) : parameter for RBF kernel
    #         dir (kxm) array: of directions
    #Ouputs: data: (2*num_sim x m array): of all the generated EC curves
    # create the GP object.
    # m here is the dimension of Euclidean space we are working in.
    data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]) )
    for (i in 1:num_sim){
        #new_data=generate_clustered_points_with_means(data_1,n=5,eta=0)
        m=generate_gp(grid_size,data_1,length_scale = length_scale)
        complex <- MatrixtoSimplicialComplex(m)
        ec_curve <- matrix(NA,nrow = 1,ncol=0)
        for (j in 1:dim(dir)[1]){
            vertex_function <- cbind(complex$Vertices[,1] , complex$Vertices%*%c(0,dir[j,1],dir[j,2],dir[j,3]))
            curve <- compute_discrete_ec_curve(complex, vertex_function, curve_length-1, first_column_index = TRUE)
            curve <- integrate_ec_curve(curve)
            # omit the length data, for now
            ec_curve <- c(ec_curve,curve[,2])
        }
        data <- rbind(data,c(1,ec_curve))
    }
    print('On class 2')
    #Second Class
    for (i in 1:num_sim){
        #new_data=generate_clustered_points_with_means(data_2,n=5,eta=0)
        m=generate_gp(grid_size,data_2,length_scale = length_scale)
        complex <- MatrixtoSimplicialComplex(m)
        ec_curve <- matrix(NA,nrow = 1,ncol=0)
        for (j in 1:dim(dir)[1]){
            vertex_function <- cbind(complex$Vertices[,1] , complex$Vertices%*%c(0,dir[j,1],dir[j,2],dir[j,3]))
            curve <- compute_discrete_ec_curve(complex, vertex_function, curve_length-1, first_column_index = TRUE)
            curve <- integrate_ec_curve(curve)
            # omit the length data, for now
            ec_curve <- c(ec_curve,curve[,2])
        }
        data <- rbind(data,c(0,ec_curve))
    }
    return(data)
}

#### Data for Generating GRF's ####
sigmoid <- function(x){
    1/(1 + exp(-x))
}
create_data <- function(num_sim, directions_specified = TRUE, dir, num_directions = 20,curve_length = 50){
    # create the GRF object
    grid_size = 20
    xy <- expand.grid(1:grid_size, 1:grid_size)
    names(xy) <- c("x","y")
    
    #specify directions for the SECT
    if (directions_specified == TRUE){
        #check if directions are there
        if(is.null(dir)) {warning("specify the directions for the SECT")}
    } else{
        dir <- generate_equidistributed_points(num_directions)
    }
    #Generate Rotation Matrices
    rotations=list()
    for (j in 1:dim(dir)[1]){
        rotation=rotation_matrix(dir[j,],c(1,0,0))
        rotations[[j]]=rotation
    }
    g.dummy1 <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, 
                      model=vgm(psill=0.025,model="Gau",range=5), nmax=20)
    yy1 <- predict(g.dummy1, newdata=xy, nsim=num_sim)
    data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]) )
    
    
    
    for(i in 3:dim(yy1)[2]){
        m <- matrix(yy1[,i],nrow = grid_size, byrow = TRUE)
        complex <- MatrixtoSimplicialComplex(m)
        ec_curve <- matrix(NA,nrow = 1,ncol=0)
        for (j in 1:dim(dir)[1]){
            vertex_function <- cbind(complex$Vertices[,1] , complex$Vertices%*%c(0,dir[j,1],dir[j,2],dir[j,3]))
            curve <- compute_discrete_ec_curve(complex, vertex_function, curve_length-1, first_column_index = TRUE)
            curve <- integrate_ec_curve(curve)
            
            # omit the length data, for now
            ec_curve <- c(ec_curve,curve[,2])
        }
        data <- rbind(data,c(1,ec_curve))
    }
    
    ## second class of random fields, given by gaussian variograms with range 6.
    g.dummy2 <- gstat(formula=z~1, locations=~x+y, dummy=T, beta=1, 
                      model=vgm(psill=0.025,model="Gau",range=6), nmax=20)
    yy2 <- predict(g.dummy2, newdata=xy, nsim=num_sim)
    
    for(i in 3:dim(yy2)[2]){
        m <- matrix(yy2[,i],nrow = grid_size, byrow = TRUE)
        complex <- MatrixtoSimplicialComplex(m)
        ec_curve <- matrix(NA,nrow = 1,ncol=0)
        for (j in 1:dim(dir)[1]){
            vertex_function <- cbind(complex$Vertices[,1] , complex$Vertices%*%c(0,dir[j,1],dir[j,2],dir[j,3]))
            curve <- compute_discrete_ec_curve(complex, vertex_function, curve_length-1, first_column_index = TRUE)
            curve <- integrate_ec_curve(curve)
            
            # omit the length data, for now
            ec_curve <- c(ec_curve,curve[,2])
        }
        data <- rbind(data,c(0,ec_curve))
    }
    data
}

#### Generating KDE data ####
generate_kde=function(grid_size=50,data,num_samples=6,num_noise=6,sd=0.1,sd_noise=0.05){
  n1=runif(num_noise,min=-1,max=1)
  n2=runif(num_noise,min=-1,max=1)
  x1=rtruncnorm(num_samples-floor(num_samples/2),a=-1,b=1,mean=data[1],sd=sd)
  y1=rtruncnorm(num_samples-floor(num_samples/2),a=-1,b=1,mean=data[2],sd=sd)
  x2=rtruncnorm(num_samples-2,a=-1,b=1,mean=0,sd=sd_noise)
  y2=rtruncnorm(num_samples-2,a=-1,b=1,mean=0,sd=sd_noise)
  points=cbind(x1,y1)
  points_2=cbind(x2,y2)
  noise=cbind(n1,n2)
  entry=rbind(points,noise)
  entry=rbind(points,points_2)
  posterior_matrix=kde(entry,gridsize=grid_size,xmin=c(-1,-1,-20),xmax=c(1,1,20))$estimate
  return(posterior_matrix)
}
create_data_kde <- function(num_sim, data_1,data_2, grid_size=20,dir,curve_length = 50,num_samples=6,num_noise=6,sd1=0.1,sd2=0.1){
  #Inputs:  num_sim: observations for each class
  #         data_1 (nxm array): the 'causal' points for class 1
  #         data_2 (nxm array) the 'causal' points for class 2
  #         grid_size (int) : number of training points on [0,1]
  #         length_scale (int) : parameter for RBF kernel
  #         dir (kxm) array: of directions
  #Ouputs: data: (2*num_sim x m array): of all the generated EC curves
  # create the GP object.
  # m here is the dimension of Euclidean space we are working in.
  sd_noise=(sd1+sd2)/2
  data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]) )
  for (i in 1:num_sim){
    #new_data=generate_clustered_points_with_means(data_1,n=5,eta=0)
    m=generate_kde(grid_size=grid_size,data=data_1,num_samples = num_samples/2,num_noise=num_noise,sd=sd1,sd_noise = sd_noise)
    complex <- MatrixtoSimplicialComplex(m)
    complex$Vertices[,2]=complex$Vertices[,2]/(grid_length/2) #?
    complex$Vertices[,3]=complex$Vertices[,3]/(grid_length/2)
    #complex$Vertices[,4]=complex$Vertices[,4]*-1
    ec_curve <- matrix(NA,nrow = 1,ncol=0)
    for (j in 1:dim(dir)[1]){
      vertex_function <- cbind(complex$Vertices[,1] , complex$Vertices%*%c(0,dir[j,1],dir[j,2],dir[j,3]))
      curve <- compute_discrete_ec_curve(complex, vertex_function, curve_length-1, first_column_index = TRUE)
      curve <- integrate_ec_curve(curve)
      # omit the length data, for now
      ec_curve <- c(ec_curve,curve[,2])
    }
    data <- rbind(data,c(1,ec_curve))
  }
  print('On class 2')
  #Second Class
  for (i in 1:num_sim){
    #new_data=generate_clustered_points_with_means(data_2,n=5,eta=0)
    m=generate_kde(grid_size,data=data_2,num_samples = num_samples*2,num_noise=num_noise*2,sd=sd2,sd_noise = sd_noise)
    complex <- MatrixtoSimplicialComplex(m)
    complex$Vertices[,2]=complex$Vertices[,2]/(grid_length/2) #?
    complex$Vertices[,3]=complex$Vertices[,3]/(grid_length/2)
    #complex$Vertices[,4]=complex$Vertices[,4]*-1
    ec_curve <- matrix(NA,nrow = 1,ncol=0)
    for (j in 1:dim(dir)[1]){
      vertex_function <- cbind(complex$Vertices[,1] , complex$Vertices%*%c(0,dir[j,1],dir[j,2],dir[j,3]))
      curve <- compute_discrete_ec_curve(complex, vertex_function, curve_length-1, first_column_index = TRUE)
      curve <- integrate_ec_curve(curve)
      # omit the length data, for now
      ec_curve <- c(ec_curve,curve[,2])
    }
    data <- rbind(data,c(0,ec_curve))
  }
  return(data)
}

#### Poisson Process and RBF Interpolation ####
#Intensity Function

intensity_func1=function(x,y){
  return(exp(5*x))
}
intensity_func2=function(x,y){
  return(exp(5*y))
}

#Random height
generate_samples_from_pp=function(pp,mean=1,sd=0.1){
  xy=cbind(pp$x,pp$y)
  xy=cbind(xy,rnorm(dim(xy)[1],mean,sd))
  return(xy)
}
difference=function(x,y){
  sumxy=sum((x-y)^2)
  return(sqrt(sumxy))
}
#RBF Gaussian
rbf_gauss=function(x,y,eta=5){
  return(exp(-(difference(x,y))))
}
#RBF inverse quadratic
inv_quad=function(x,y,eta=5){
  return(1/(sqrt(1+(difference(x,y)/eta)^2)))
}
#RBF computation
compute_rbf_model=function(data,eta=5,func){
  x=data[,1:2]
  y=data[,3]
  samples=dim(x)[1]
  kernel_matrix=matrix(NA,ncol=samples,nrow=samples)
  for (i in 1:samples){
    kernel_matrix[i,]=apply(X=x,FUN=func,MARGIN=1,y=x[i,])
  }
  lambdas=solve(kernel_matrix)%*%y
  rbf_final=list(lambdas=lambdas,x=x,values=y,eta=eta,func=func)
}
#RBF prediction
rbf_predict=function(rbf_model,grid){
  lambdas=rbf_model$lambdas
  x=rbf_model$x
  func=rbf_model$func
  num_preds=dim(grid)[1]
  num_centers=dim(x)[1]
  predictions=rep(0,num_preds)
  for (i in 1:num_centers){
    distances=apply(X = grid,FUN = func,MARGIN=1,y=x[i,])
    distances=distances*lambdas[i]
    predictions=predictions+distances
  }
  return(predictions)
}
rbf_on_grid=function(grid_size=25,func=rbf_gauss,data,eta=5){
  rbf_model=compute_rbf_model(data=data,eta=eta,func=func)
  grids=seq(-1,1,length=grid_size)
  grid=expand.grid(grids,grids)
  predictions=rbf_predict(rbf_model = rbf_model,grid=grid)
  predictions=matrix(predictions,nrow=grid_size)
  return(predictions)
}

generate_poisson_complex=function(grid_size=25,func=rbf_gauss,eta=5,intensity_func,poisson_parameter=4){
  #Poisson Process
  pp1=rpoispp(intensity_func,poisson_parameter,win=owin(c(-1,1),c(-1,1)))
  if (pp1$n<2){
    return(generate_poisson_complex(grid_size=grid_size,func=func,eta=eta,intensity_func=intensity_func,poisson_parameter=poisson_parameter))
  }
  samples1=generate_samples_from_pp(pp1,sd=0.25)
  #### Interpolation ####
  predictions1=rbf_on_grid(grid_size=grid_size,func=rbf_gauss,data=samples1,eta=eta)
  complex1=MatrixtoSimplicialComplexTriangular(predictions1,grid_length=grid_size)
  complex=list(complex=complex1,base_points=samples1)
  return(complex)
}

create_data_poisson=function(num_sim=25,dir,curve_length=10,poisson_parameter=4,grid_size=25,func=rbf_gauss,eta=5){
  data <- matrix(NA,nrow=0,ncol = 1+curve_length*( dim(dir)[1]) )
  for (i in 1:num_sim){
    complex=generate_poisson_complex(grid_size = grid_size,func = func,eta = eta,intensity_func = intensity_func1,poisson_parameter = poisson_parameter)
    complex=complex[[1]]
    ec_curve <- matrix(NA,nrow = 1,ncol=0)
    for (j in 1:dim(dir)[1]){
      vertex_function <- complex$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])
      curve <- compute_discrete_ec_curve(complex, vertex_function, curve_length-1, first_column_index = FALSE)
      curve <- integrate_ec_curve(curve)
      # omit the length data, for now
      ec_curve <- c(ec_curve,curve[,2])
    }
    data <- rbind(data,c(1,ec_curve))
  }
  print('On class 2')
  #Second Class
  for (i in 1:num_sim){
    complex=generate_poisson_complex(grid_size = grid_size,func = func,eta = eta,intensity_func = intensity_func2,poisson_parameter = poisson_parameter)
    complex=complex[[1]]
    ec_curve <- matrix(NA,nrow = 1,ncol=0)
    for (j in 1:dim(dir)[1]){
      vertex_function <- complex$Vertices%*%c(dir[j,1],dir[j,2],dir[j,3])
      curve <- compute_discrete_ec_curve(complex, vertex_function, curve_length-1, first_column_index = FALSE)
      curve <- integrate_ec_curve(curve)
      # omit the length data, for now
      ec_curve <- c(ec_curve,curve[,2])
    }
    data <- rbind(data,c(0,ec_curve))
  }
  return(data)
}
#### Helper Functions for Plotting and SECT ###

#Convert the matrix to a simplicial complex (4 vertices per face)
MatrixtoSimplicialComplex <- function(matrix){
    vertices <- matrix(NA,nrow = 0, ncol = 4)
    edges <- matrix(NA, nrow = 0, ncol = 2)
    faces <- matrix(NA,nrow = 0, ncol = 4)
    length_x <- dim(matrix)[1]
    length_y <- dim(matrix)[2]
    
    for(i in 1:(length_x - 1)){
        for(j in 1:(length_y - 1)){
            #The last three columns are the coordinates of the vertices
            #indices of the vertices of this pixel:
            vertex1 <- (i-1)*(dim(matrix)[2]) + j
            vertex2 <- (i - 1)*(dim(matrix)[2]) + j+1
            vertex3 <- (i)*(dim(matrix)[2]) + j
            vertex4 <- (i)*(dim(matrix)[2]) + j+1
            
            # add the vertices to complex
            vertices <- rbind(vertices, c( vertex1, i-length_x/2, j-length_y/2, matrix[i,j]) )
            vertices <- rbind(vertices, c( vertex2, i-length_x/2, j+1-length_y/2, matrix[i,j+1]))
            vertices <- rbind(vertices, c( vertex3, i+1-length_x/2, j-length_y/2, matrix[i+1,j] ))
            vertices <- rbind(vertices, c( vertex4, i+1-length_x/2, j+1-length_y/2, matrix[i+1,j+1] ))
            
            # add the edges to complex
            edges <- rbind(edges, c(vertex1,vertex2))
            edges <- rbind(edges, c(vertex2,vertex4))
            edges <- rbind(edges, c(vertex1,vertex3))
            edges <- rbind(edges, c(vertex3,vertex4))
            # add the faces to complex
            faces <- rbind(faces, c(vertex1,vertex2,vertex3,vertex4))
            
        }
    }
    
    vertices <- unique(vertices)
    vertices = vertices[order(vertices[,1]),]
    edges <- unique(edges)
    faces <- unique(faces)
    
    simplicial_complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
}
#Convert the matrix to a triangular simplicial complex (3 vertices per face)
MatrixtoSimplicialComplexTriangular <- function(matrix,grid_length){
    vertices <- matrix(NA,nrow = 0, ncol = 4)
    edges <- matrix(NA, nrow = 0, ncol = 2)
    faces <- matrix(NA,nrow = 0, ncol = 3)
    length_x <- dim(matrix)[1]
    length_y <- dim(matrix)[2]
    
    for(i in 1:(length_x - 1)){
        for(j in 1:(length_y - 1)){
            #The last three columns are the coordinates of the vertices
            #indices of the vertices of this pixel:
            vertex1 <- (i-1)*(dim(matrix)[2]) + j
            vertex2 <- (i - 1)*(dim(matrix)[2]) + j+1
            vertex3 <- (i)*(dim(matrix)[2]) + j
            vertex4 <- (i)*(dim(matrix)[2]) + j+1
            
            # add the vertices to complex
            vertices <- rbind(vertices, c( vertex1, i-length_x/2, j-length_y/2, matrix[i,j]) )
            vertices <- rbind(vertices, c( vertex2, i-length_x/2, j+1-length_y/2, matrix[i,j+1]))
            vertices <- rbind(vertices, c( vertex3, i+1-length_x/2, j-length_y/2, matrix[i+1,j] ))
            vertices <- rbind(vertices, c( vertex4, i+1-length_x/2, j+1-length_y/2, matrix[i+1,j+1] ))
            
            # add the edges to complex
            edges <- rbind(edges, c(vertex1,vertex2))
            edges <- rbind(edges, c(vertex2,vertex4))
            edges <- rbind(edges, c(vertex1,vertex3))
            edges <- rbind(edges, c(vertex3,vertex4))
            #edges<-rbind(edges,c(vertex2,vertex3))
            # add the faces to complex
            #faces <- rbind(faces, c(vertex1,vertex2,vertex3,vertex4))
            faces <- rbind(faces, c(vertex1,vertex2,vertex3))
            faces <- rbind(faces, c(vertex2,vertex3,vertex4))
        }
    }
    
    vertices <- unique(vertices)
    vertices = vertices[order(vertices[,1]),]
    edges <- unique(edges)
    faces <- unique(faces)
    complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
    complex$Vertices[,2]=complex$Vertices[,2]/(grid_length/2)
    complex$Vertices[,3]=complex$Vertices[,3]/(grid_length/2)
    complex$Vertices=complex$Vertices[,2:4]
    vertex=complex$Vertices
    faces=complex$Faces
    ind=rep(1,dim(vertex)[1])
    vertex=cbind(vertex,ind)
    mesh=tmesh3d(t(vertex),indices=t(faces))
    vertices=as.matrix(t(mesh$vb)[,1:3])
    faces=as.matrix(t(mesh$it))
    edges=vcgGetEdge(mesh)
    edges=as.matrix(edges[,1:2])
    complex <- list(Vertices = vertices, Edges = edges, Faces = faces)
    return(complex)
}
reshape_complex=function(complex,grid_length){
  complex$Vertices[,2]=complex$Vertices[,2]/(grid_length/2)
  complex$Vertices[,3]=complex$Vertices[,3]/(grid_length/2)
  complex$Vertices=complex$Vertices[,2:4]
  return(complex)
}
convert_complex=function(complex){
  #posterior_complex2$Vertices[,4]=posterior_complex2$Vertices[,4]/(-1)
  vertex=complex$Vertices
  faces=complex$Faces
  ind=rep(1,dim(vertex)[1])
  vertex=cbind(vertex,ind)
  mesh=tmesh3d(t(vertex),indices=t(faces))
  return(mesh)
}

