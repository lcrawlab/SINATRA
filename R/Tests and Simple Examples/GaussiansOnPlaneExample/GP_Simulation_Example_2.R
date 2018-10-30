source('../Final Code for Functions (Really)/Critical_Points_and_plot.R')
source('../Final Code for Functions (Really)/GP_Generation.R')
source('../Final Code for Functions (Really)/SECT_final.R')
source('Generate_fields.R')

# how about start with flat plane with added normal densities add several points
### Generating Two Classes of GPs ###
# specify the training points that are to be fixed.
noncausal_points = matrix(c(0,0,-0.5,
                            0.1,0.1,0.5,
                            0.2,0.2,-0.5,
                            0.3,0.3,0.5,
                            0.4,0.4,-0.5,
                            0.5,0.5,0.5,
                            0.6,0.6,-0.5,
                            0.7,0.7,0.5,
                            0.8,0.8,-0.5), nrow = 9, byrow = TRUE)

#specify the class specific points.
class_one_causal_points = matrix(c(0.75,0.25,-0.5,
                                   0.95,0.1,0.5,
                                   0.6,0.4,-0.5),nrow = 3, byrow = TRUE)

class_two_causal_points = matrix(c(0.25,0.75,-0.5,
                                   0.1,0.95,0.5,
                                   0.4,0.6,-0.5),nrow = 3, byrow = TRUE)

class_one_points = rbind(noncausal_points, class_one_causal_points) # recenter the points
class_two_points = rbind(noncausal_points, class_two_causal_points)

#generate the desired directions
#directions = generate_equidistributed_points_hemisphere(25)
directions = matrix(c(0,0,1, 1/sqrt(2),0,1/sqrt(2), 0,1/sqrt(2),1/sqrt(2)), ncol = 3, byrow = TRUE)

grid_length = 50
data = create_data_gp(num_sim = 50,grid_size = grid_length, dir = directions)

#### Visualize some GPs from each class ####

posterior_matrix1=generate_gp(grid_length,class_one_points,length_scale = 5)
posterior_complex1 <- MatrixtoSimplicialComplexTriangular(posterior_matrix1)
posterior_complex1$Vertices[,2]=posterior_complex1$Vertices[,2]/(grid_length/2) #?
posterior_complex1$Vertices[,3]=posterior_complex1$Vertices[,3]/(grid_length/2)

#create the simplicial complex
try_vertex1=posterior_complex1$Vertices[,2:4]
try_faces1=posterior_complex1$Faces
ind1=rep(1,dim(try_vertex1)[1])
try_vertex1=cbind(try_vertex1,ind1)
grf_mesh_posterior1=tmesh3d(t(try_vertex1),indices=t(try_faces1))
plot3d(grf_mesh_posterior1,col = alpha.col('grey',alpha=0.1))
aspect3d(1,1,1)

filename='GP_1'
#Write to one OFF file
vcgOffWrite(grf_mesh_posterior1,filename=filename)

#Repeat the same thing with the second class of GP
posterior_matrix2=generate_gp(grid_length,class_two_points,length_scale = 5)
posterior_complex2 <- MatrixtoSimplicialComplexTriangular(posterior_matrix2)
posterior_complex2$Vertices[,2]=posterior_complex2$Vertices[,2]/(grid_length/2)
posterior_complex2$Vertices[,3]=posterior_complex2$Vertices[,3]/(grid_length/2)
try_vertex2=posterior_complex2$Vertices[,2:4]
try_faces2=posterior_complex2$Faces
ind2=rep(1,dim(try_vertex2)[1])
try_vertex2=cbind(try_vertex2,ind2)
grf_mesh_posterior2=tmesh3d(t(try_vertex2),indices=t(try_faces2))
plot3d(grf_mesh_posterior,colors=colors)
aspect3d(1,1,1)

filename='GP_2'
#Write to one OFF file
vcgOffWrite(grf_mesh_posterior2,filename=filename)


### Find Critical Points in examples from each class ###
grf1=process_off_file_v2('GP_1.off')
phi=0.1
len=50-1
#critical_point_coordinates1=as.data.frame(crit_points_coordinates_general(dir[1,],grf1,len,phi))
critical_points_total1=critical_points_multiple_directions(dir=directions,off=grf1,steps=len,phi=phi)
critical_point_coordinates1=get_all_critical_points(critical_points_total1)
#Get the critical points that correspond to the selected region
#Now for GRF2
grf2=process_off_file_v2('GP_2.off')
#critical_point_coordinates2=as.data.frame(crit_points_coordinates_general(dir[1,],grf2,len,phi))
critical_points_total2=critical_points_multiple_directions(dir=directions,off=grf2,steps=len,phi=phi)
critical_point_coordinates2=get_all_critical_points(critical_points_total2)




### Feature Selection and Reconstruction ###
# We use Rate. Fit a Gaussian process classifier to this data
gp <- gausspr(x = train_data,y = train_labels, kernel = 'rbfdot', kpar = list(sigma = 0.015))
label_probabilities <- predict(gp, test_data, type = 'probabilities')

predicted_labels <- ifelse(label_probabilities > 0, 1,0)

cm = table(predicted_labels, test_labels)
cm

correct = (cm[1,1] + cm[2,2])/sum(cm)
correct

n <- dim(gp_data)[1]
X <- gp_data[,-1]
h <- 0.01 #median(dist(X)) #need a better choice of this; how does bandwidth affect kernel choice?

f <- rep(0,n)
Kn <- GaussKernel(t(X),1/(2*h^2))
diag(Kn)=1

# do Newton IRLS procedure to use Laplace / Gaussian approximation.
for(k in 1:600){
  W <- diag(as.vector(sigmoid(f)*(1-sigmoid(f))))
  B <- diag(x = 1,n) + sqrt(W) %*% Kn %*% sqrt(W)
  #Kinda show
  L <- chol(B)
  b <- W%*%f + (gp_data[,1]+1)/2 - sigmoid(f)
  #Kinda slow
  #L has components that 
  #Try to not use all the solves
  a <- b - solve(sqrt(W)%*%t(L),solve(L,sqrt(W)%*%Kn%*%b))
  f <- Kn%*%a
}
v = solve(L, sqrt(W)%*%Kn)
# generate samples from approximate posterior
#fhat samples may not be right dimensions
fhat.samples = rmvnorm(1e4,f, Kn - t(v)%*% v)


# use RATE:
cores = cores=detectCores()
#Dimension of Res should be 50 samples
res = RATEv2(X=X,f.draws=fhat.samples,prop.var = 1,snp.nms = colnames(X),cores = cores)

#What we want
rates=res$RATE
fhat.samples=0
B=0
W=0
L=0
f=0
Kn=0
want=rates>1/length(rates)
numeric_want=as.numeric(want)
want_indices=which(1==numeric_want)

#### Finding Critical Points ####
#Rate's Critical Points
selected_crit_points1_rate=as.matrix(critical_points_interval_multiple_directions(grf1,
                                                                              critical_points_total1,want_indices,len,direction = directions))
selected_crit_points2_rate=as.matrix(critical_points_interval_multiple_directions(grf2,
                                                                              critical_points_total2,want_indices,len,direction = dirrections))
#### Plotting Time####
#Plotting RATE's Critical Points 
colors1=rep('white',dim(grf_mesh_posterior1$vb)[2])
colors2=rep('white',dim(grf_mesh_posterior2$vb)[2])

vertex_transpose1=t(grf_mesh_posterior1$vb)[,1:3]
vertex_transpose2=t(grf_mesh_posterior2$vb)[,1:3]


plot_image_crit_points(as.matrix(selected_crit_points1_rate),vertex_transpose1,n=9,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior1)
find_n_closest_vertices(pretty_new_means,vertices=vertex_transpose1,n=9,color = 'blue')
#find_n_closest_vertices(noise,vertices=vertex_transpose1,n=9,color='green')
open3d()
aspect3d(1,1,1)
plot_image_crit_points(as.matrix(selected_crit_points2_rate),vertex_transpose2,n=9,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior2)
find_n_closest_vertices(pretty_new_means2,vertices=vertex_transpose2,n=5,color = 'blue')
#find_n_closest_vertices(noise,vertices=vertex_transpose2,n=5,color='green')
close3d()










