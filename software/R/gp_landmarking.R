#### Rough Implementation of GP Landmarking ####
#setwd('/Users/brucewang/Dropbox (DataPlusMath)/Data + Experiments Tim Sudijono/')
library(Rvcg)
library(Matrix)
library(FNN)
library(pracma)
# This is basically an R implentation of getting Gaussian Process Landmarks

#### Compute Curvatures ####


compute_curvature = function(vertices, k){
  points = t(vertices)
  n = knnx.index(points,points,k+1)
  n = n[,-1]
  p = repmat(points,k,1) - points[as.vector(n),]
  p = array(p,dim = c(dim(points)[1],k,3))
  nv = dim(points)[1]
  
  c = matrix(0,nrow = nv, ncol = 6)
  
  c[,1] = rowSums(p[,,1] * p[,,1])
  c[,2] = rowSums(p[,,1] * p[,,2])
  c[,3] = rowSums(p[,,1] * p[,,3])
  c[,4] = rowSums(p[,,2] * p[,,2])
  c[,5] = rowSums(p[,,2] * p[,,3])
  c[,6] = rowSums(p[,,3] * p[,,3])
  
  c = c/k
  curvature = rep(0,dim(points)[1])
  for (i in 1:nv){
    cmat = matrix(c(c[i,1],c[i,2],c[i,3],c[i,2],c[i,4],c[i,5],c[i,3],c[i,5],c[i,6]), nrow = 3)
    v = eigen(cmat)$vectors
    d = eigen(cmat)$values
    
    lambda = min(d)
    k = which(d == lambda)
    curvature[i] = lambda / sum(d)
  }
  return(curvature)
}

find_landmarks = function(mesh, num_landmark = 10){
  
  curv = compute_curvature(mesh$vb[-4,],10)
  
  nv = nverts(mesh)
  nf = nfaces(mesh)
  
  I = cbind(matrix(mesh$it[1,], nrow = 1),matrix(mesh$it[2,],nrow = 1), matrix(mesh$it[3,],nrow =1))
  J = matrix(rep(1:nf, 3),nrow = 1)
  f2v = sparseMatrix(i = J,j = I,x = 1,dims = c(nf,nv))
  
  face_area = vcgArea(mesh, TRUE)$pertriangle
  
  vert_area = (face_area%*% f2v)/3
  
  ep = 0.5
  #lambda = vert_area * (ep * abs(gauss_curv)/sum(abs(gauss_curv))) + (1-ep) * abs(avg_curv)/sum(abs(avg_curv))
  lambda = vert_area * curv/sum(curv)
  
  verts = mesh$vb[1:3,]
  edges = vcgGetEdge(mesh)
  edges1 = edges[,1]
  edges2 = edges[,2]
  bandwidth = mean(sqrt(colSums((verts[,edges1] - verts[,edges2])^2)))/5
  BNN = min(500, nv)
  #Computing the 500 KNN, rest is considered to be 0
  idx = knnx.index(t(verts),t(verts),BNN+1)
  dist = knnx.dist(t(verts),t(verts),BNN+1)
  #Distance Matrix
  full_phi = sparseMatrix(i = rep(1:nv,BNN+1),j = idx,x = as.vector(exp(-dist^2/bandwidth)),dims = c(nv,nv))
  full_phi = (full_phi+t(full_phi))/2
  
  #print('Constructing Full Kernel')
  
  full_mat_prod = full_phi %*% sparseMatrix(i = (1:nv), j = 1:nv,x = lambda,dims = c(nv,nv)) %*% full_phi
  
  kernel_trace = diag(full_mat_prod)
  
  
  landmark_index = rep(0,num_landmark)
  inv_kn = matrix(0,ncol = num_landmark, nrow = num_landmark)
  for (k in 1:num_landmark){
    if (k == 1){
      ptuq = kernel_trace
    }
    else{
      if (k == 2){
        inv_kn[1:(k-1),1:(k-1)] = 1/full_mat_prod[landmark_index[1],landmark_index[1]]
        ptuq = kernel_trace - colSums(t(full_mat_prod[,landmark_index[1:(k-1)]]) 
                                      * (inv_kn[1:(k-1),1:(k-1)]) %*% t(full_mat_prod[landmark_index[1:(k-1)],]))
      }
      else{
        p = full_mat_prod[landmark_index[1:(k-2)], landmark_index[k-1]]
        mu = 1/ (full_mat_prod[landmark_index[k-1],landmark_index[k-1]] - t(p) %*% inv_kn[1:(k-2),1:(k-2)]%*%p)
        inv_kn[1:(k-2),1:(k-1)] = inv_kn[1:(k-2),1:(k-2)] %*% cbind(diag(k-2) + mu[1,1] * (matrix(p,nrow = length(p)) %*% t(p)) * inv_kn[1:(k-2),1:(k-2)], (-1*mu[1,1]) * p)
        inv_kn[(k-1),1:(k-1)] = c(t(inv_kn[1:(k-2),k-1]),mu)
        productEntity = inv_kn[1:(k-1),1:(k-1)]%*%full_mat_prod[landmark_index[1:(k-1)],]
        ptuq = t(kernel_trace - colSums(t(full_mat_prod[,landmark_index[1:(k-1)]]) * productEntity))
      }
    }
    max_index = which(ptuq == max(ptuq))
    landmark_index[k] = max_index
  }
  return(landmark_index)
}

plot_selected_landmarks = function(mesh, selected_vertices, num_landmarks = 20){
  vertices = t(mesh$vb[-4,])
  selected_vertices = c(0,selected_vertices)
  landmarks = c(0,find_landmarks(mesh, num_landmarks))
  lmk1 = intersect(selected_vertices,landmarks)
  lmk2 = landmarks[-(which(landmarks %in% lmk1))]
  found_landmarks = find_landmarks_from_vertices(vertices, selected_vertices, landmarks)
  lmk1 = found_landmarks$found_landmarks
  lmk2 = found_landmarks$lost_landmarks
 # print(lmk1)
 # print(lmk2)
  vertices[,3] = vertices[,3] - 0.01
  rgl.points(vertices[lmk1,] ,col = 'orange', size = 7)
  rgl.points(vertices[lmk2,] ,col = 'black', size = 7)
}

find_landmarks_from_vertices = function(vertices, selected_vertices, landmarks, n = 10){
  found_landmarks = c()
  lost_landmarks = c()
  closest_points = knnx.index(vertices, vertices, k = n)
  for (i in 1:length(landmarks)){
    landmark = landmarks[i]
    if (landmark %in% closest_points[selected_vertices,]){
      found_landmarks = c(found_landmarks,landmark)
    }
    else{
      lost_landmarks = c(lost_landmarks,landmark)
    }
  }
  return(list(found_landmarks = found_landmarks, lost_landmarks = lost_landmarks))
}

