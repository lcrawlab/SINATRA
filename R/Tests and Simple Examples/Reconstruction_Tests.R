source('Final Code for Functions (Really)/Critical_Points_and_plot.R')
source('Final Code for Functions (Really)/Feature_Selection.R')
source('Final Code for Functions (Really)/SECT_final.R')

##### Let's create a simple test case: simplicial complex representing a sphere #####
# Read OFF file
sphere <- process_off_file_v3('Topological_Reconstruction/sphere_off_file.off')
colors=rep('white',dim(sphere$Vertices)[1])
#mr_sphere=vcgImport('Topological_Reconstruction/sphere_off_file.off')
# visualize
sphere_obj <- visualize_OFF(sphere)
#plot3d(mr_sphere,col=colors)

# choose a single direction, compute the SECT + critical points
directions <- matrix(c(0,0,1),nrow = 1, byrow = TRUE)

all_critical_points_sphere = critical_points_multiple_directions(dir=directions,off=sphere,steps=50,phi=0.1)
critical_point_coordinates_sphere = get_all_critical_points(all_critical_points_sphere)
plot_image_crit_points(points=critical_point_coordinates_sphere, vertices = sphere$Vertices, n=5, off_file=sphere_obj)
#plot_image_crit_points(points=critical_point_coordinates_sphere, vertices = sphere$Vertices, n=5, off_file=mr_sphere)


# Choose several directions this time.

##### Let's create a simple test case: simplicial complex representing a sphere #####
# Read OFF file
dragon <- process_off_file_v3('Topological_Reconstruction/dragon_off_file.off')
mr_dragon=vcgImport('Topological_Reconstruction/dragon_off_file.off')
colors=rep('white',dim(dragon$Vertices)[1])
plot3d(mr_dragon,col=colors)
# visualize
dragon_obj <- visualize_OFF(dragon)

# choose a single direction, compute the SECT + critical points
directions <- matrix(c(0,0,1),nrow = 1, byrow = TRUE)
directions <- matrix(c(1,0,0),nrow = 1, byrow = TRUE)
directions <- matrix(c(0,0,1,1,0,0,0,1,0),ncol = 3, byrow = TRUE)
directions <- matrix(c(2,0.1,0),ncol = 3, byrow = TRUE)
directions <- matrix(c(1/sqrt(3),1/sqrt(3),1/sqrt(3)),nrow = 1, byrow = TRUE)


# compute the critical points
all_critical_points_dragon = critical_points_multiple_directions(dir=directions,off=dragon,steps=5,phi=0.001)
critical_point_coordinates_dragon = get_all_critical_points(all_critical_points_dragon)
critical_point_coordinates_dragon
#plot_image_crit_points(points=critical_point_coordinates_dragon, vertices = dragon$Vertices, n=5, off_file=dragon_obj,
#                   axes=TRUE,labels = TRUE,plot_points = TRUE, hull = FALSE)
#plot_image_crit_points(points=critical_point_coordinates_dragon, vertices = dragon$Vertices, n=5, off_file=mr_dragon,axes=TRUE,labels=TRUE,plot_points = TRUE)

points = matrix(critical_point_coordinates_dragon,ncol = 3,byrow = FALSE)
vertices = dragon$Vertices
off_file=dragon_obj

number_crit_points=dim(points)[1]


# draw a small point at that place
colors=rep('white',dim(off_file$vb)[2])
plot3d(off_file,xlab='',ylab='',zlab='',type=c('shade'),col=colors,alpha=0.5,axes=FALSE,labels=FALSE)

# Show tick marks
axis3d('x', pos=c( NA, 0, 0 ), col = "red")
axis3d('y', pos=c( 0, NA, 0 ), col = "green")
axis3d('z', pos=c( 0, 0, NA ), col = "black")
for (i in 1:number_crit_points){
  dummy_point=points[i,]
  dummy_point=matrix(dummy_point,ncol=3)
  closest_points_index=knnx.index(vertices,dummy_point,k=1,algo = 'kd_tree') # find the closest points to the off vertices
  closest_points=matrix(vertices[closest_points_index,],ncol=3,byrow = TRUE)
  rgl::points3d(closest_points[,1],closest_points[,2],closest_points[,3],col='blue',alpha=0.5,size = 10)
  rgl::points3d(dummy_point[,1],dummy_point[,2],dummy_point[,3],col='red',alpha=0.5,size = 10)
  #rgl.points(dummy_point,color=color,alpha=1,size=10)
}

### Molars ####
molar <- process_off_file_v3('Data/HDM/by_diet/Follivore/MCZ-5070_M1080.off')
colors=rep('white',dim(molar$Vertices)[1])
# visualize
molar_obj <- visualize_OFF(molar)

# choose a single direction, compute the SECT + critical points
directions <- matrix(c(0,0,1),nrow = 1, byrow = TRUE)
directions <- matrix(c(1,0,0),nrow = 1, byrow = TRUE)
directions <- matrix(c(0,0,1,1,0,0,0,1,0),ncol = 3, byrow = TRUE)
directions <- matrix(c(2,0.1,0),ncol = 3, byrow = TRUE)
directions <- matrix(c(1/sqrt(3),1/sqrt(3),1/sqrt(3)),nrow = 1, byrow = TRUE)


# Compute the EC curves
vertex_function=molar$Vertices%*%c(0,0,1)
ec_curve=compute_discrete_ec_curve(molar,vertex_function,len)
ec_curve



# compute the critical points
all_critical_points_molar = critical_points_multiple_directions(dir=directions,off=molar,steps=50,phi=0.001)
critical_point_coordinates_molar = get_all_critical_points(all_critical_points_molar)
critical_point_coordinates_molar
#plot_image_crit_points(points=critical_point_coordinates_dragon, vertices = dragon$Vertices, n=5, off_file=dragon_obj,
#                   axes=TRUE,labels = TRUE,plot_points = TRUE, hull = FALSE)
#plot_image_crit_points(points=critical_point_coordinates_dragon, vertices = dragon$Vertices, n=5, off_file=mr_dragon,axes=TRUE,labels=TRUE,plot_points = TRUE)


points = matrix(critical_point_coordinates_molar,ncol = 3,byrow = FALSE)
vertices = molar$Vertices
off_file=molar_obj

number_crit_points=dim(points)[1]


# draw a small point at that place
colors=rep('white',dim(off_file$vb)[2])
plot3d(off_file,xlab='',ylab='',zlab='',type=c('shade'),col=colors,alpha=0.5,axes=FALSE,labels=FALSE)

# Show tick marks
axis3d('x', pos=c( NA, 0, 0 ), col = "red")
axis3d('y', pos=c( 0, NA, 0 ), col = "green")
axis3d('z', pos=c( 0, 0, NA ), col = "black")
for (i in 1:number_crit_points){
  dummy_point=points[i,]
  dummy_point=matrix(dummy_point,ncol=3)
  closest_points_index=knnx.index(vertices,dummy_point,k=1,algo = 'kd_tree') # find the closest points to the off vertices
  closest_points=matrix(vertices[closest_points_index,],ncol=3,byrow = TRUE)
  rgl::points3d(closest_points[,1],closest_points[,2],closest_points[,3],col='blue',alpha=0.5,size = 10)
  rgl::points3d(dummy_point[,1],dummy_point[,2],dummy_point[,3],col='red',alpha=0.5,size = 10)
  #rgl.points(dummy_point,color=color,alpha=1,size=10)
}


### Helper Functions ###

visualize_OFF <- function(simp_complex){
  vertices=simp_complex$Vertices[,1:3]
  faces=simp_complex$Faces
  indices = rep(1,dim(vertices)[1])
  vertices = cbind(vertices,indices)

  off_obj=tmesh3d(t(vertices),indices=t(faces))
  plot3d(off_obj,col = alpha.col('grey',alpha=0.1))
  off_obj
}

############# Poisson Process Example #############

source('Final Code for Functions (Really)/Critical_Points_and_plot.R')
source('Final Code for Functions (Really)/SECT_final.R')
source('Final Code for Functions (Really)/GP_Generation.R')
library(spatstat)
####Generate Poisson Processs on Grid####


pp1=rpoispp(intensity_func1,4,win=owin(c(-1,1),c(-1,1)))
pp2=rpoispp(intensity_func2,4,win=owin(c(-1,1),c(-1,1)))
samples1=generate_samples_from_pp(pp1,sd=0.25)
samples2=generate_samples_from_pp(pp2,sd=0.25)


#### Interpolation ####
eta=1
grid_size=50
predictions1=rbf_on_grid(grid_size=grid_size,func=rbf_gauss,data=samples1,eta=eta)
complex1=MatrixtoSimplicialComplexTriangular(predictions1,grid_length=grid_size)
newcomplex1=convert_complex(complex1)
plot_image_crit_points(samples1,vertices = complex1$Vertices,n=2,axes=TRUE,off_file=newcomplex1,hull=FALSE)
predictions2=rbf_on_grid(grid_size=grid_size,func=rbf_gauss,data=samples2,eta=eta)
complex2=MatrixtoSimplicialComplexTriangular(predictions2,grid_length = grid_size)
newcomplex2=convert_complex(complex2)
open3d()
plot_image_crit_points(samples2,vertices = complex2$Vertices,n=2,axes=TRUE,off_file=newcomplex2,hull=FALSE)

#alpha_complex=alphaComplexFiltration(complex1$Vertices)
#### Critical Point Computation ####


dir=generate_equidistributed_points_hemisphere(50)
phi=0.0001
len=10

critical_points_total1=critical_points_multiple_directions(dir=dir,off=complex1,steps=len,phi=phi)
#vertex_function=complex1$Vertices%*%dir[4,]
#ec_curve=compute_discrete_ec_curve(complex1,vertex_function,len)
critical_point_coordinates1=get_all_critical_points(critical_points_total1)
#critical_point_coordinates1=critical_point_coordinates1[-23,]
#critical_point_coordinates1[4,]=c(0,0,1.0334)
#shite_points= crit_points_coordinates(list(c(1,0,0)),off=complex1,steps=50,eta=0.1,delta=0.1)
print(critical_points_total1)
#open3d()
plot_image_crit_points(critical_point_coordinates1,vertices=complex1$Vertices,off_file=newcomplex1,hull=FALSE,axes=TRUE,point_size=5)
plot_points(points=samples1,point_color = 'red',point_size = 7)
plot_image_crit_points(critical_point_coordinates1,vertices=complex1$Vertices,n=10,off_file=newcomplex1,axes=TRUE,point_size=5)

critical_points_total2=critical_points_multiple_directions(dir=dir,off=complex2,steps=len,phi=phi)
#vertex_function=complex1$Vertices%*%dir[4,]
#ec_curve=compute_discrete_ec_curve(complex1,vertex_function,len)
critical_point_coordinates2=get_all_critical_points(critical_points_total2)
print(critical_points_total2)
#open3d()
plot_image_crit_points(n=10,critical_point_coordinates2,vertices=complex2$Vertices,off_file=newcomplex2,axes=TRUE,point_size=5)
plot_points(points=samples2,point_color = 'red',point_size = 7)



