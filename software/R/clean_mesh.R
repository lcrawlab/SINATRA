#### Aligning and Scaling the meshes ####
#library(Morpho)
#library(rgl)
#library(Rvcg)

#' Scale & Translate Meshes
#'
#' @export
#'
#' @import Rvcg
#' @import Morpho
#' @import rgl
#'
#' @description Scales the meshes to have unit surface area and translates the meshes to the origin.
#'
#' @param input_dir The input directory of the files to scale & translate.
#' @param output_dir The output diirectory to put the files post scaling and translating.
#' @return None; files written in output path.
clean_files=function(input_dir,output_dir){
    files=list.files(path = dir,full.names = TRUE)
    filenames=list.files(path=dir,full.names = FALSE)
    num_files=length(files)
    for (i in 1:num_files){
        print(files[i])
        file=vcgImport(files[i])
        filename=paste(out_path,filenames[i],sep='')
        area=vcgArea(file)
        file2=scalemesh(file,1/sqrt(area),center='mean')
        centroid <- colMeans(vert2points(file2))
        file3=translate3d(file2,-centroid[1],-centroid[2],-centroid[3])
        print(filename)
        vcgOffWrite(file3,filename = filename)
    }
}



