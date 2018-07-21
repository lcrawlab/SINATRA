library(Rcpp)
library(RcppArmadillo)
library(mvtnorm)
library(kernlab)
library(doParallel)
setwd('/Users/Bruce/Desktop/Research_Code/TDA')
#source('SECT_3D_tim.R')
#source('Finding_Critical_Points.R')
source('SECT_Final.R')
source('Critical_Points_and_plot.R')
source('GP_Generation.R')
sourceCpp("BAKRGibbs.cpp")
source('RATEv2.R')
library(truncnorm)
library(varbvs)

#### Generate the Data ####
eff_seed <- sample(1:2^15, 1)
print(eff_seed)
set.seed(eff_seed)

n1=rtruncnorm(5,a=-1,b=1,mean=0,sd=1)
n2=rtruncnorm(5,a=-1,b=1,mean=0,sd=1)
n3=rtruncnorm(5,a=-0.1,b=0.1,mean=0,sd=0.1)
x1=rtruncnorm(5,a=-0.5,b=1,mean=0.35,sd=0.2)
y1=rtruncnorm(5,a=-0.5,b=1,mean=0.2,sd=0.25)
z1=rtruncnorm(5,a=-1,b=1,mean=0.24,sd=0.2)
x2=rtruncnorm(5,a=-1,b=0.5,mean=-0.3,sd=0.15)
y2=rtruncnorm(5,a=-1,b=0.5,mean=-0.2,sd=0.25)
z2=rtruncnorm(5,a=-1,b=1,mean=-0.24,sd=0.1)
noise=as.matrix(cbind(n1,n2,n3))
pretty_new_means=as.matrix(cbind(x1,y1,z1))
data_1=as.matrix(rbind(pretty_new_means,noise))
#pretty_new_means=as.matrix(rbind(pretty_new_means,noise))
#keep_index=sample(10,5,replace=FALSE)
#same_indices=pretty_new_means[keep_index,]
pretty_new_means2=as.matrix(cbind(x2,y2,z2))
data_2=as.matrix(rbind(pretty_new_means2,noise))
#pretty_new_means2=as.matrix(rbind(pretty_new_means2,noise))
grid_len=20
#dir = matrix(c(0,0,1),nrow = 1)
dir = matrix(c(0,0,1, 1/sqrt(2),0,1/sqrt(2), 0,1/sqrt(2),1/sqrt(2)), ncol = 3, byrow = TRUE)
gp_data=create_data_gp(num_sim = 100, data_1 = data_1,data_2 = data_2, grid_size=20,length_scale=5,dir= dir)

ind <- sample(1:dim(gp_data)[1],50)
train_data <- gp_data[-ind,-1]
train_labels <- gp_data[-ind,1]
test_data <- gp_data[ind,-1]
test_labels <- gp_data[ind,1]

#### Generate Two GP's (from each class) #### 
post_matrix1=generate_gp(grid_len,data_1,length_scale = 5)
posterior_complex1 <- MatrixtoSimplicialComplexTriangular(post_matrix1)
posterior_complex1$Vertices[,2]=posterior_complex1$Vertices[,2]/(grid_len/2)
posterior_complex1$Vertices[,3]=posterior_complex1$Vertices[,3]/(grid_len/2)
try_vertex1=posterior_complex1$Vertices[,2:4]
try_faces1=posterior_complex1$Faces
ind1=rep(1,dim(try_vertex1)[1])
try_vertex1=cbind(try_vertex1,ind1)
grf_mesh_posterior1=tmesh3d(t(try_vertex1),indices=t(try_faces1))
#plot3d(grf_mesh_posterior,colors=colors)
#aspect3d(1,1,1)
filename='GP_1'
#Write to one OFF file
vcgOffWrite(grf_mesh_posterior1,filename=filename)
data_2=generate_clustered_points_with_means(pretty_new_means2,n=5,eta=0)

post_matrix2=generate_gp(grid_len,data_2,length_scale = 5)
posterior_complex2 <- MatrixtoSimplicialComplexTriangular(post_matrix2)
posterior_complex2$Vertices[,2]=posterior_complex2$Vertices[,2]/(grid_len/2)
posterior_complex2$Vertices[,3]=posterior_complex2$Vertices[,3]/(grid_len/2)
try_vertex2=posterior_complex2$Vertices[,2:4]
try_faces2=posterior_complex2$Faces
ind2=rep(1,dim(try_vertex2)[1])
try_vertex2=cbind(try_vertex2,ind2)
grf_mesh_posterior2=tmesh3d(t(try_vertex2),indices=t(try_faces2))
#plot3d(grf_mesh_posterior,colors=colors)
#aspect3d(1,1,1)
filename='GP_2'
#Write to one OFF file
vcgOffWrite(grf_mesh_posterior2,filename=filename)

####Finding Critical Points ####
grf1=process_off_file_v3('GP_1.off')
phi=0.1
len=50-1
#critical_point_coordinates1=as.data.frame(crit_points_coordinates_general(dir[1,],grf1,len,phi))
critical_points_total1=critical_points_multiple_directions(dir=dir,off=grf1,steps=len,phi=phi)
critical_point_coordinates1=get_all_critical_points(critical_points_total1)
#Get the critical points that correspond to the selected region
#Now for GRF2
grf2=process_off_file_v3('GP_2.off')
#critical_point_coordinates2=as.data.frame(crit_points_coordinates_general(dir[1,],grf2,len,phi))
critical_points_total2=critical_points_multiple_directions(dir=dir,off=grf2,steps=len,phi=phi)
critical_point_coordinates2=get_all_critical_points(critical_points_total2)



#### Feature Selection ####
### RATE ###
#Fit a Gaussian process classifier to this data
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

#RATE
f <- rep(0,n)
Kn <- GaussKernel(t(X),1/(2*h^2))
diag(Kn)=1

#Change stopping conditions for convergence

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


library(matrixcalc)
#is.positive.semi.definite(q_hat)
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


#### Lasso ####

library(glmnet)
regression_model=cv.glmnet(gp_data[,-1], gp_data[,1], alpha = 0.95,intercept = FALSE)
label_probabilities <- predict(regression_model, test_data)

predicted_labels <- ifelse(label_probabilities > 0, 1,0)
cm = table(predicted_labels, test_labels)
cm
correct = (cm[1,1] + cm[2,2])/sum(cm)
correct
#Extract the coefficients
tmp_coeffs = as.matrix(coef(regression_model,s=regression_model$lambda.1se))
regression_coeff=as.vector(tmp_coeffs)
regression_coeff=regression_coeff[2:length(regression_coeff)]
want_indices_lasso=which(0<regression_coeff)

#### Bayesian Variable Selection ####
fit=varbvs(X = gp_data[,-1],Z = NULL, y= gp_data[,1],family='gaussian')
label_probabilities <- predict(fit, test_data)

predicted_labels <- ifelse(label_probabilities > 0, 1,0)
cm = table(predicted_labels, test_labels)
cm
correct = (cm[1,1] + cm[2,2])/sum(cm)
correct
bayesian_probs=fit$pip
want_indices_bayesian=which(0.05<bayesian_probs)

#### Finding Critical Points ####
#Rate's Critical Points
selected_crit_points1_rate=as.matrix(critical_points_interval_multiple_directions(grf1,critical_points_total1,want_indices,len,direction = dir))
selected_crit_points2_rate=as.matrix(critical_points_interval_multiple_directions(grf2,critical_points_total2,want_indices,len,direction = dir))
#PPA's Critical Points
selected_crit_points1_bayesian=as.matrix(critical_points_interval_multiple_directions(grf1,critical_points_total1,want_indices_bayesian,len,direction = dir))
selected_crit_points2_bayesian=as.matrix(critical_points_interval_multiple_directions(grf2,critical_points_total2,want_indices_bayesian,len,direction = dir))
#Lasso Critical Points
selected_crit_points1_lasso=as.matrix(critical_points_interval_multiple_directions(grf1,critical_points_total1,want_indices_lasso,len,direction = dir))
selected_crit_points2_lasso=as.matrix(critical_points_interval_multiple_directions(grf2,critical_points_total2,want_indices_lasso,len,direction = dir))

num_crit_points_grf1=c(dim(selected_crit_points1_rate)[1],dim(selected_crit_points1_bayesian)[1],dim(selected_crit_points1_lasso)[1],dim(critical_point_coordinates1)[1])
num_crit_points_grf2=c(dim(selected_crit_points2_rate)[1],dim(selected_crit_points2_bayesian)[1],dim(selected_crit_points2_lasso)[1],dim(critical_point_coordinates2)[1])
number_critical_points=rbind(num_crit_points_grf1,num_crit_points_grf2)
colnames(number_critical_points)=c('Rate','PPA','Lasso','Total Critical Points')
number_critical_points

#### Plotting Time####
#Plotting RATE's Critical Points 
colors1=rep('white',dim(grf_mesh_posterior1$vb)[2])
colors2=rep('white',dim(grf_mesh_posterior2$vb)[2])

vertex_transpose1=t(grf_mesh_posterior1$vb)[,1:3]
vertex_transpose2=t(grf_mesh_posterior2$vb)[,1:3]


plot_image_crit_points(as.matrix(selected_crit_points1_rate),vertex_transpose1,n=9,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior1)
find_n_closest_vertices(pretty_new_means,vertices=vertex_transpose1,n=9,color = 'blue')
find_n_closest_vertices(noise,vertices=vertex_transpose1,n=9,color='green')
open3d()
aspect3d(1,1,1)
plot_image_crit_points(as.matrix(selected_crit_points2_rate),vertex_transpose2,n=9,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior2)
find_n_closest_vertices(pretty_new_means2,vertices=vertex_transpose2,n=5,color = 'blue')
find_n_closest_vertices(noise,vertices=vertex_transpose2,n=5,color='green')
close3d()

#Plotting PPA's Critical Points 

plot_image_crit_points(as.matrix(selected_crit_points1_bayesian),vertex_transpose1,n=9,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior1)
find_n_closest_vertices(pretty_new_means,vertices=vertex_transpose1,n=9,color = 'blue')
find_n_closest_vertices(noise,vertices=vertex_transpose1,n=9,color='green')
open3d()
aspect3d(1,1,1)
plot_image_crit_points(as.matrix(selected_crit_points2_bayesian),vertex_transpose2,n=9,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior2)
find_n_closest_vertices(pretty_new_means2,vertices=vertex_transpose2,n=9,color = 'blue')
find_n_closest_vertices(noise,vertices=vertex_transpose2,n=9,color='green')
close3d()

#Plotting Lasso's Critical Points 

plot_image_crit_points(as.matrix(selected_crit_points1_lasso[,-1]),vertex_transpose1,n=9,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior1)
find_n_closest_vertices(pretty_new_means,vertices=vertex_transpose1,n=9,color = 'blue')
find_n_closest_vertices(noise,vertices=vertex_transpose1,n=9,color='green')
open3d()
aspect3d(1,1,1)
plot_image_crit_points(as.matrix(selected_crit_points2_lasso),vertex_transpose2,n=9,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior2)
find_n_closest_vertices(pretty_new_means2,vertices=vertex_transpose2,n=9,color = 'blue')
find_n_closest_vertices(noise,vertices=vertex_transpose2,n=9,color='green')
close3d()

#Plotting all critical points

plot_image_crit_points(as.matrix(critical_point_coordinates1),vertex_transpose1,n=5,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior1)
find_n_closest_vertices(pretty_new_means,vertices=vertex_transpose1,n=9,color = 'blue')
find_n_closest_vertices(noise,vertices=vertex_transpose1,n=9,color='green')
open3d()
aspect3d(1,1,1)
plot_image_crit_points(as.matrix(critical_point_coordinates2),vertex_transpose2,n=5,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior2)
find_n_closest_vertices(pretty_new_means2,vertices=vertex_transpose2,n=9,color = 'blue')
find_n_closest_vertices(noise,vertices=vertex_transpose2,n=9,color='green')
close3d()
#plot_image_crit_points(as.matrix(critical_point_coordinates2),vertex_transpose1,n=20,alpha_one=0.75,crit_color='red',off_file = grf_mesh2)


plot_image_crit_points(as.matrix(critical_points_total1[[3]][,-1]),vertex_transpose1,n=5,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior1)
find_n_closest_vertices(pretty_new_means,vertices=vertex_transpose1,n=9,color = 'blue')
find_n_closest_vertices(noise,vertices=vertex_transpose1,n=9,color='green')
open3d()
aspect3d(1,1,1)
plot_image_crit_points(as.matrix(critical_points_total2[[3]][,-1]),vertex_transpose2,n=5,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior2)
find_n_closest_vertices(pretty_new_means2,vertices=vertex_transpose2,n=9,color = 'blue')
find_n_closest_vertices(noise,vertices=vertex_transpose2,n=9,color='green')

#### Test to get Critical Points ####
#x1=rnorm(10,0.1,0.4)
#y1=rnorm(10,0,0.2)
#z1=rnorm(10,0,0.4)
n1=rtruncnorm(5,a=-1,b=1,mean=0,sd=1)
n2=rtruncnorm(5,a=-1,b=1,mean=0,sd=15)
n3=rtruncnorm(5,a=-0.1,b=0.1,mean=0,sd=0.1)
x1=rtruncnorm(5,a=-0.5,b=1,mean=0.35,sd=0.2)
y1=rtruncnorm(5,a=-0.5,b=1,mean=0.2,sd=0.25)
z1=rtruncnorm(5,a=-1,b=1,mean=0.24,sd=0.2)
x2=rtruncnorm(5,a=-1,b=0.5,mean=-0.3,sd=0.15)
y2=rtruncnorm(5,a=-1,b=0.5,mean=-0.2,sd=0.25)
z2=rtruncnorm(5,a=-1,b=1,mean=-0.24,sd=0.1)
noise=as.matrix(cbind(n1,n2,n3))
pretty_new_means=as.matrix(cbind(x1,y1,z1))
pretty_new_means1=as.matrix(rbind(pretty_new_means,noise))
#pretty_new_means=as.matrix(rbind(pretty_new_means,noise))
#keep_index=sample(10,5,replace=FALSE)
#same_indices=pretty_new_means[keep_index,]
pretty_new_means2=as.matrix(cbind(x2,y2,z2))
pretty_new_means2.5=as.matrix(rbind(pretty_new_means2,noise))
#pretty_new_means2=as.matrix(rbind(pretty_new_means2,noise))
grid_len=20
dir = matrix(c(0,0,1),nrow = 1)
#noisy_data=generate_clustered_points_with_means(noise,n=5,eta=0)
#clustered_1=generate_clustered_points_with_means(pretty_new_means,n=5,eta=0)
#data_1=as.matrix(rbind(noisy_data,clustered_1))
post_matrix1=generate_gp(grid_len,pretty_new_means1,length_scale = 5)
#post_matrix1=generate_gp(grid_len,pretty_new_means,length_scale = 5)
posterior_complex1 <- MatrixtoSimplicialComplexTriangular(post_matrix1)
posterior_complex1$Vertices[,2]=posterior_complex1$Vertices[,2]/(grid_len/2)
posterior_complex1$Vertices[,3]=posterior_complex1$Vertices[,3]/(grid_len/2)
try_vertex1=posterior_complex1$Vertices[,2:4]
try_faces1=posterior_complex1$Faces
ind1=rep(1,dim(try_vertex1)[1])
try_vertex1=cbind(try_vertex1,ind1)
grf_mesh_posterior1=tmesh3d(t(try_vertex1),indices=t(try_faces1))
#plot3d(grf_mesh_posterior,colors=colors)
#aspect3d(1,1,1)
filename='GP_1'
#Write to one OFF file
vcgOffWrite(grf_mesh_posterior1,filename=filename)
#clustered_2=generate_clustered_points_with_means(pretty_new_means2,n=5,eta=0)
#data_2=as.matrix(rbind(noisy_data,clustered_2))
post_matrix2=generate_gp(grid_len,pretty_new_means2.5,length_scale = 5)
#post_matrix2=generate_gp(grid_len,pretty_new_means2,length_scale = 5)
posterior_complex2 <- MatrixtoSimplicialComplexTriangular(post_matrix2)
posterior_complex2$Vertices[,2]=posterior_complex2$Vertices[,2]/(grid_len/2)
posterior_complex2$Vertices[,3]=posterior_complex2$Vertices[,3]/(grid_len/2)
try_vertex2=posterior_complex2$Vertices[,2:4]
try_faces2=posterior_complex2$Faces
ind2=rep(1,dim(try_vertex2)[1])
try_vertex2=cbind(try_vertex2,ind2)
grf_mesh_posterior2=tmesh3d(t(try_vertex2),indices=t(try_faces2))
#plot3d(grf_mesh_posterior,colors=colors)
#aspect3d(1,1,1)
filename='GP_2'
#Write to one OFF file
vcgOffWrite(grf_mesh_posterior2,filename=filename)
grf1=process_off_file_v3('GP_1.off')
direcs=list()
direcs[[1]]=c(0,0,1)
eta=0.01
delta=0.01
len=50-1
critical_point_coordinates1=as.data.frame(find_critical_points_coordinates_modified(direcs,grf1,len,eta,delta))
#Get the critical points that correspond to the selected region
ec_curve=extract_ec_curve(direction=c(0,0,1),grf1,len)
#Now for GRF2
grf2=process_off_file_v3('GP_2.off')
critical_point_coordinates2=as.data.frame(find_critical_points_coordinates_modified(direcs,grf2,len,eta,delta))
ec_curve2=extract_ec_curve(direction=c(0,0,1),grf2,len)

vertex_transpose1=t(grf_mesh_posterior1$vb)[,1:3]
vertex_transpose2=t(grf_mesh_posterior2$vb)[,1:3]
plot_image_crit_points(as.matrix(critical_point_coordinates1),vertex_transpose1,n=5,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior1)
find_n_closest_vertices(pretty_new_means,vertices=vertex_transpose1,n=9,color = 'blue')
find_n_closest_vertices(noise,vertices=vertex_transpose1,n=9,color='green')
open3d()
aspect3d(1,1,1)
plot_image_crit_points(as.matrix(critical_point_coordinates2),vertex_transpose2,n=5,alpha_one=0.75,crit_color='red',off_file = grf_mesh_posterior2)
find_n_closest_vertices(pretty_new_means2,vertices=vertex_transpose2,n=9,color = 'blue')
find_n_closest_vertices(noise,vertices=vertex_transpose2,n=9,color='green')

#### Finding EC Curves ###

complex <- MatrixtoSimplicialComplex(post_matrix1)
vertex_function <- cbind(complex$Vertices[,1] , complex$Vertices%*%c(0,dir[1,1],dir[1,2],dir[1,3]))
curve <- compute_discrete_ec_curve(complex, vertex_function, 49, first_column_index = TRUE)
curve <- integrate_ec_curve(curve)
# omit the length data, for now
rotation=rotation_matrix(c(1,1,1),c(0,0,1))
complex2=complex
new_vert=cbind(complex2$Vertices[,1],complex2$Vertices[,2:4]%*%rotation)
complex2$Vertices=new_vert
vertex_function2=cbind(complex2$Vertices[,1] , complex2$Vertices%*%c(0,dir[1,1],dir[1,2],dir[1,3]))
curve2 <- compute_discrete_ec_curve(complex2, vertex_function2, 49, first_column_index = TRUE)
curve2 <- integrate_ec_curve(curve2)
new_dir=c(1,1,1)
vertex_function3 <- cbind(complex$Vertices[,1] , complex$Vertices%*%c(0,new_dir[1],new_dir[2],new_dir[3]))
curve3 <- compute_discrete_ec_curve(complex, vertex_function3, 49, first_column_index = TRUE)
curve3 <- integrate_ec_curve(curve3)

