###### 
# This script tries out the idea of having an evidence parameter for each critical point based on how many projections are selected by the critical point algorithm.
# The end goal is heatmap on the shape of which points are the most important
######
source('R/Reconstruction_Functions.R')
source('R/SECT_Functions.R')
source('R/Simulation_Functions.R')
source('R/Feature_Selection_Functions.R')

### First, generate two classes of Poisson Point Processes
grid_size=25; func=rbf_gauss; poisson_parameter=4

dir=generate_equidistributed_points(100)
phi=0.001
#Curve Length
len=5
nsim=35
eta=0.1
#Generate Data
data=create_data_poisson(num_sim=nsim,curve_length=len,poisson_parameter = poisson_parameter,grid_size=grid_size,dir=dir,eta=eta)
#Feature Selection with rate, select the features only if they're above 1/p, where p is the number of features
rate_features=find_rate_variables(data,radius=1,bandwidth=0.01,weights = FALSE)
#rate_features0=find_rate_variables(data,radius=0,bandwidth=0.01)

#Generate example poisson processes
poisson1=generate_poisson_complex(grid_size = grid_size,func = func,eta = eta,intensity_func = intensity_func1,poisson_parameter = poisson_parameter)
complex1=poisson1[[1]]
newcomplex1=convert_complex(complex1)
#Class 2
poisson2=generate_poisson_complex(grid_size = grid_size,func = func,eta = eta,intensity_func = intensity_func2,poisson_parameter = poisson_parameter)
complex2=poisson2[[1]]
newcomplex2=convert_complex(complex2)
