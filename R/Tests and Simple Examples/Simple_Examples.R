#Make sure to install the packages needed! There is a lot. 
source('Final Code for Functions (Really)/Reconstruction_Functions.R')
source('Final Code for Functions (Really)/SECT_final.R')
#Examples of Critical Point finding for simple Shapes

#The direction (You can modify this to play around)
dir <- matrix(c(1,0,0),ncol=3, byrow = TRUE)
#The curve length (you don't need to worry about it too much)
curve_length=50
#Angle of rotation (also not too important)
phi=0.1
####Tetrahedron####
tetra=process_off_file_v3('Data/Simple_Shapes/tetra.off')
tetra_plot=vcgImport('Data/Simple_Shapes/tetra.off')
#Find the critical Points
#### Example ####
tetra
a=generate_perturbed_directions(direction = c(1,0,0),phi=0.1)
a
find_critical_points_curve(direction=c(1,0,0),tetra,50,phi=0.1)
find_critical_point_coordinates(direction = c(1,0,0),off=tetra,steps=50,phi=0.1,rejection_sampling = FALSE)
#### End example 
critical_points_tetra=get_all_critical_points(find_critical_points_multiple_directions(dir =dir,off=tetra,steps = curve_length,phi = phi,rejection_sampling = FALSE))
#Plot the critical points for tetrahedron
plot_complex_with_points(points = critical_points_tetra,hull = FALSE,crit_color = 'blue',off_file = tetra_plot)

#### Cube ####
cube=process_off_file_v3('Data/Simple_Shapes/cube.off')
cube_plot=vcgImport('Data/Simple_Shapes/cube.off')
#Find the critical Points
critical_points_cube=get_all_critical_points(find_critical_points_multiple_directions(dir =dir,off=cube,steps = curve_length,phi = phi))
#Plot the critical points for tetrahedron
plot_complex_with_points(points = critical_points_cube,hull = FALSE,off_file = cube_plot)


#### Docedahedron ####
dodec=process_off_file_v3('Data/Simple_Shapes/dodec.off')
dodec_plot=vcgImport('Data/Simple_Shapes/dodec.off')
#Find the critical Points
critical_points_dodec=get_all_critical_points(find_critical_points_multiple_directions(dir =dir,off=dodec,steps = curve_length,phi = phi))
#Plot the critical points for tetrahedron
plot_complex_with_points(points = critical_points_dodec,hull = FALSE,off_file = dodec_plot)


#### Bean ####
dir <- matrix(c(0,0,1),ncol=3, byrow = TRUE)
bean=process_off_file_v3('Data/Simple_Shapes/bean_off_file.off')
bean_plot=vcgImport('Data/Simple_Shapes/bean_off_file.off')
#Find the critical Points
critical_points_bean=get_all_critical_points(find_critical_points_multiple_directions(dir =dir,off=bean,steps = curve_length,phi = phi,rejection_sampling = FALSE))
#Plot the critical points for tetrahedron
plot_complex_with_points(points = critical_points_bean,hull = FALSE,off_file = bean_plot)


#### Mushroom ####
mushroom=process_off_file_v3('Data/Simple_Shapes/mushroom.off')
mushroom_plot=vcgImport('Data/Simple_Shapes/mushroom.off')
#Find the critical Points
critical_points_mushroom=get_all_critical_points(find_critical_points_multiple_directions(dir =dir,off=mushroom,steps = curve_length,phi = phi))
#Plot the critical points for tetrahedron
plot_complex_with_points(points = critical_points_mushroom,hull = FALSE,off_file = mushroom_plot)
