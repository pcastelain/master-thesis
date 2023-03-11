##################################################################
##                     Dependence Inference                     ##
##################################################################

#WARNING: You should only run this script once the marginal parameters have already been estimated (cf script n�1)

#Load marginal estimates
load("../Tmp/Storm_Marginals.RData")

#Load functions to generate larger SPDE domain and project locations
source("Functions/utils_domain_and_projections.R")
#Load functions used in the inference process
source("Functions/utils_inference.R")
#Load functions to plot maps, obtain resulting covariance and draw random samples from the r-Pareto process
source("Functions/utils_reporting.R")


#A few functions need to be sent to all workers in parallel optimization. They are packaged together in a package called "nonStatInf"
#Install last version of the packaged functions
library(devtools)
setwd("nonStatInf/")
install()
library(nonStatInf)
packageDescription("nonStatInf")
setwd("../")


library(EnvStats)
library(ggmap)
library(parallel)
library(patchwork)

##################################################################
##                             Maps                             ##
##################################################################

zoom_factor <- 1.1
bbox <- c(min(loc.sub[,1])-(max(loc.sub[,1])-min(loc.sub[,1]))*(zoom_factor-1),min(loc.sub[,2])-(max(loc.sub[,2])-min(loc.sub[,2]))*(zoom_factor-1),max(loc.sub[,1])+(max(loc.sub[,1])-min(loc.sub[,1]))*(zoom_factor-1),max(loc.sub[,2])+(max(loc.sub[,2])-min(loc.sub[,2]))*(zoom_factor-1))

#Here we use a boundary box not to use the zoom functionality
map <- get_map(location = bbox, source="stamen", maptype = "toner" , color="bw") #A map mainly blakc in the background ! Usefull to draw heatmap
map_lite <- get_map(location = bbox, source="stamen", maptype = "toner-lite") #A map mainly white in background. Used to draw vector fields


#################################################################
##                         SPDE domain                         ##
#################################################################


#Generate a larger domain with 30% of buffer on each side.
big_domain_loc <- generate_larger_domain(loc.sub,0.3)

#Projection of the big domain
big_domain_proj <- loc_to_proj(big_domain_loc,fixed_lat = mean(loc.sub[,2]))

#indices of our domain locations within the big domain
big_domain_rest_indices <- get_region_indices(big_domain_loc,loc.sub)

list_init <- init_anis(length(unique(big_domain_loc[,1])),length(unique(big_domain_loc[,2])),as.data.frame(big_domain_proj),restriction_indices = big_domain_rest_indices,land_sea = land_sea)

save(big_domain_loc,big_domain_proj,big_domain_rest_indices,list_init,file="../Tmp/list_init.RData")
load("../Tmp/list_init.RData")


#################################################################
##                   Adjust weight functions                   ##
#################################################################

weight_coeff <- 0.1

keep_looping <- TRUE

while(keep_looping){
  rescaled_storm_weight <- rep(NA, length(storm_footprint[,1]))
  for(i in c(1:length(rescaled_storm_weight))){
    rescaled_storm_weight[i] <- (1 - exp((1 - c(rescaled_storm_footprint[i,]%*%r_func_vec))*weight_coeff))
  }
  if(quantile(rescaled_storm_weight,0.3) < 0.9){ #we want 30% < 0.9 not more
    weight_coeff <- weight_coeff + 0.1
  }else{
    print(paste0("Weight_Coeff: ",weight_coeff))
    keep_looping <- FALSE
  }
}


weighting_function <- function(x, weight_coeff,r_func_vec){
  x * (1 - exp((1 - c(x%*%r_func_vec))*weight_coeff))
}
weighting_function_derivative <- function(x, weight_coeff, r_func_vec){
  (1 - exp((1 - c(x%*%r_func_vec))*weight_coeff)) + (x*r_func_vec) * exp((1 - c(x%*%r_func_vec))*weight_coeff)
}

#Number of location in composite gradient score and number of composite
subset_estimation_size <- 25
nb_composite <- 100 #incresing the number of composite likelihood increases the robustness of the results. You might use 1000

list_index <- lapply(1:nb_composite, function(i){sample(c(1:length(storm_footprint[1,])),subset_estimation_size)})


list_storm_footprint <- lapply(1:length(storm_footprint[,1]), function(i){rescaled_storm_footprint[i,]})##################here I used the rescaled ones !



#Demonstration with a polynomial model of degree 1 ! Intialization values can be found in appendix of the thesis report for other models if needed
initial_values <- c(0.000944068, 52.969863403, 6.167713939, 2.927578822, 0.006520219, -0.005668306, 0.001715515, -0.003897556)

param_scale <- c(0.001,10,rep(1,2),rep(10e-3,4)) #rescaling the parameters can significantly help the algorithm converge

library(marqLevAlg)
estimate_MLA<- marqLevAlg(b = c(initial_values)/param_scale, #The argument m=8 could be used if we have no initial values
                           fn = objective_function_gradient_score_anis,
                           excesses = list_storm_footprint,
                           anis_list = list_init,
                           weighting_function = weighting_function,
                           weighting_function_derivative = weighting_function_derivative,
                           weight_coeff = weight_coeff,
                           r_func_vec = r_func_vec,
                           complete_list_anis = complete_list_anis,
                           compos_like_index_list = list_index,
                            param_scale = param_scale, #If a scaling is used it needs to be provided to the score function !
                           anis_mode = 1, #This is the anisotropy model used, 1 for polynomial and 2 for interpolation cf. implementation of "complete_list_anis" and the thesis report
                           print.info = TRUE,
                           nproc = detectCores(),#Here this can be changed 
                           maxiter=10,
                           .packages = c("nonStatInf","akima","Matrix"),
                           file="../Debug/debug_mla.txt")

estimate_MLA$b <- estimate_MLA$b * param_scale
save(estimate_MLA,file="../Tmp/saved_estimate_MLA.RData")
load("../Tmp/saved_estimate_MLA.RData")

#Covariance matrix associated with the fitted model
covarianceMatrix <- get_cov_mat(complete_list_anis(list_init, parameter = estimate_MLA$b, anis_mode = 1))
#List of anisotropic functions and SPDE parameters associated with the fitted model 
list_complete <- complete_list_anis(list_init, parameter = estimate_MLA$b,anis_mode = 1)




##################################################################
##                         Extremogramm                         ##
##################################################################


# Compute all possible pairwise combinations
pairs_combination_indexes <- expand.grid(1:length(loc.sub[,1]),1:length(loc.sub[,1])) #First proposition mais on va faire plus leger car sinon les points sont dupliqués
#pairs_combination_indexes <- t(combn(1:length(loc.sub[,1]),2)) #Comme tout ce que on calcule est symetrique on s'en fou ! Ca diminue par 2 la taille si on compare au précédent !


# Create vector of empirical extremogram
empirical_extremogram <- rep(NA, length(pairs_combination_indexes[,1]))

quantileLevel_storm <- 0.52 #0.625


quantile_storm <- rep(NA, length(storm_footprint[1,]))

for(i in 1:length(storm_footprint[1,])){
  quantile_storm[i] <- quantile(storm_footprint[,i], quantileLevel_storm, na.rm = TRUE)
}


# Compute empirical extremogram for each pair of locations
for(k in 1:length(empirical_extremogram)){
  empirical_extremogram[k] <- computeExtremalCoeff(pairs_combination_indexes[k,1], pairs_combination_indexes[k,2], storm_footprint, quantile_storm)
  #print(paste0(k,"/",length(empirical_extremogram)))
}


index_to_plot <- sample(length(empirical_extremogram),5000)#sample(length(empirical_extremogram),5000)

projected_coordinates <- cbind(diag(list_init$D_i)[list_init$restriction_indices],diag(list_init$D_j)[list_init$restriction_indices])#loc_to_proj(loc.sub,scaling = TRUE,scale = attr(big_domain_proj,"scaled:scale"),center = attr(big_domain_proj,"scaled:center"))

# Compute vector of distance between each pairs
distances_between_pairs <- sqrt((projected_coordinates[pairs_combination_indexes[,1],1] - projected_coordinates[pairs_combination_indexes[,2],1])^2 +
                                  (projected_coordinates[pairs_combination_indexes[,1],2] - projected_coordinates[pairs_combination_indexes[,2],2])^2)
difference_between_pairs_x <- (projected_coordinates[pairs_combination_indexes[,1],1] - projected_coordinates[pairs_combination_indexes[,2],1])
difference_between_pairs_y <- (projected_coordinates[pairs_combination_indexes[,1],2] - projected_coordinates[pairs_combination_indexes[,2],2])


# Compute orientation of each pairs
matrixOfAngles <- atan((projected_coordinates[pairs_combination_indexes[,1],1] - projected_coordinates[pairs_combination_indexes[,2],1]) /
                         (projected_coordinates[pairs_combination_indexes[,1],2] - projected_coordinates[pairs_combination_indexes[,2],2])) / pi * 180
matrixOfAngles[is.nan(matrixOfAngles)] <- 0


# Store results in data frame for easier plotting
estimates <- data.frame(extremogram = empirical_extremogram[!is.nan(empirical_extremogram)], dist = distances_between_pairs[!is.nan(empirical_extremogram)], angle = matrixOfAngles[!is.nan(empirical_extremogram)],
                        x = difference_between_pairs_x[!is.nan(empirical_extremogram)], y = difference_between_pairs_y[!is.nan(empirical_extremogram)])

# Plot empirical extremogram as function of distance
p_Extremogram_Dist <- ggplot(data = estimates[index_to_plot,], aes(x = dist, y = extremogram, color = angle)) + geom_point(alpha = 0.5) +
  guides(color = guide_colorbar(barheight = 10, barwidth = 1)) + scale_color_gradientn(colors = rev(rainbow(7))) + ylim(0,1)
print(p_Extremogram_Dist)

# Plot empirical extremogram as a map
p_Extremogram_Map <- ggplot(data = estimates[index_to_plot,], aes(x = x, y = y, color = extremogram)) + geom_point(alpha = 0.5,size=1.5) +
  guides(color = guide_colorbar(barheight = 10, barwidth = 1)) + scale_color_gradientn(colors = rev(rainbow(7)), limits = c(0,1)) + coord_fixed(ratio = 1) +
  xlab("Distance along X") + ylab("Distance along Y") #+ guide=guide_legend(title = "Extremogram")
#, limits = c(0.5,1)
print(p_Extremogram_Map)



# vario_model <- rep(0,length(pairs_combination_indexes[,1]))
# 
# for(i in 1:length(pairs_combination_indexes[,1])){
#   
#   # Compute extremogram between pairs
# 
#   s1 <- pairs_combination_indexes[i,1]
#   s2 <- pairs_combination_indexes[i,2]
# 
#   vario_model[i] <- 1/2 * (covarianceMatrix[s1,s1]+covarianceMatrix[s2,s2]) - covarianceMatrix[s1,s2]
# 
#   extremogram_model[i] <- 2 * (1 - pnorm(sqrt(vario_model[i] / 2)))
# }



# Compute theoretical extremogram from variogram
extremogram_model <- rep(0,length(pairs_combination_indexes[,1]))

##### Test pour calculer plus vite !
library(Matrix)
S1 <- sparseMatrix(j=c(1:length(pairs_combination_indexes[,1])),i=pairs_combination_indexes[,1],x=1,dims=c(length(storm_footprint[1,]),length(pairs_combination_indexes[,1])))
S2 <- sparseMatrix(j=c(1:length(pairs_combination_indexes[,1])),i=pairs_combination_indexes[,2],x=1,dims=c(length(storm_footprint[1,]),length(pairs_combination_indexes[,1])))

vario_model <- rep(0,length(pairs_combination_indexes[,1]))

split_mat <- c(round(seq(1,length(pairs_combination_indexes[,1]),by = 1000)),length(pairs_combination_indexes[,1]))
for(k in c(2:length(split_mat))){
  vario_model[split_mat[k-1]:split_mat[k]] <- diag(1/2 * (t(S1[,split_mat[k-1]:split_mat[k]]-S2[,split_mat[k-1]:split_mat[k]]) %*% Matrix(covarianceMatrix,sparse = TRUE) %*% (S1[,split_mat[k-1]:split_mat[k]]- S2[,split_mat[k-1]:split_mat[k]])))
}

extremogram_model <- 2 * (1 - pnorm(sqrt(vario_model / 2)))

model <- data.frame(extremogram = extremogram_model, dist = distances_between_pairs, angle = matrixOfAngles,x = difference_between_pairs_x[!is.nan(empirical_extremogram)], y = difference_between_pairs_y[!is.nan(empirical_extremogram)])
                    #x = distances_between_pairs[!is.nan(empirical_extremogram)] * cos(matrixOfAngles[!is.nan(empirical_extremogram)]), y = distances_between_pairs[!is.nan(empirical_extremogram)] * sin(matrixOfAngles[!is.nan(empirical_extremogram)]))


#k <- 120
#index_to_plot <- which(pairs_combination_indexes[,1] %in% k | pairs_combination_indexes[,2] %in% k )



# Plot empirical variogram cloud along with estimated model
pExtremogram <-
  ggplot(data = estimates[index_to_plot, ], aes(x = dist, y = extremogram, color = angle)) + geom_point(alpha = 0.2) + geom_smooth(se = FALSE) +
  geom_line(data = model[index_to_plot, ],
            mapping =  aes(x = dist, y = extremogram, color = angle)) +
  guides(color = guide_colorbar(barheight = 10, barwidth = 1)) + scale_color_gradientn(colors = rev(rainbow(7))) + ylim(0, 1)
print(pExtremogram)

# Plot model extremogram as a map
p_Extremogram_model_Map <- ggplot(data = model[index_to_plot,], aes(x = x, y = y, color = extremogram)) + geom_point(alpha = 0.5) +
  guides(color = guide_colorbar(barheight = 10, barwidth = 1)) + scale_color_gradientn(colors = rev(rainbow(7)), limits = c(0,1)) #, limits = c(0.5,1)
print(p_Extremogram_model_Map)


################Create extremogram maps
extremogram_model_map <- matrix(NA,ncol=length(loc.sub[,1]),nrow=length(loc.sub[,1]))
extremogram_empirical_map <- matrix(NA,ncol=length(loc.sub[,1]),nrow=length(loc.sub[,1]))
for(i in 1:length(extremogram_model)){
  
  # Compute extremogram between pairs
  
  s1 <- pairs_combination_indexes[i,1]
  s2 <- pairs_combination_indexes[i,2]
  
  extremogram_model_map[s1,s2] =  extremogram_model[i]
  extremogram_model_map[s2,s1] =  extremogram_model[i]
  
  extremogram_empirical_map[s1,s2] =  empirical_extremogram[i]
  extremogram_empirical_map[s2,s1] =  empirical_extremogram[i]
}

#################################################################
##                  Plotting extremogram maps                  ##
#################################################################
# 2 functions to easily compare modeled extremal dependence with empirical one

plot_empirical <- function(i){
  nice_plot(extremogram_empirical_map[i,],loc.sub,map,lims = c(0,1))
}
plot_model  <- function(i){
  nice_plot(extremogram_model_map[i,],loc.sub,map,lims = c(0,1))
}

plot_empirical(238) + plot_model(238)

#################################################################
##               Test best empirical pi quantile               ##
#################################################################
#The quantile with minimum squared difference between empirical and modeled extremogram can be used to set "quantileLevel_storm" at line 340
  
for(thres in seq(0.40,0.7,length.out = 20)){
  
  
  quantile_storm <- rep(NA, length(storm_footprint[1,]))
  
  for(i in 1:length(storm_footprint[1,])){
    quantile_storm[i] <- quantile(storm_footprint[,i], thres, na.rm = TRUE)
  }
  
  
  # Compute empirical extremogram for each pair of locations
  for(k in 1:length(empirical_extremogram)){
    empirical_extremogram[k] <- computeExtremalCoeff(pairs_combination_indexes[k,1], pairs_combination_indexes[k,2], storm_footprint, quantile_storm)
  }
  
  print(paste("Threshol quantile:",thres,"Square diff:",sum((extremogram_model- empirical_extremogram)^2)))
}

#################################################################
##                      Draw vector field                      ##
#################################################################

plot_vector_field <- function(list_complete,gamma=0.1,size=0.75){
  x <- list_complete$D_i#[list_complete$restriction_indices,list_complete$restriction_indices]
  y <- list_complete$D_j#[list_complete$restriction_indices,list_complete$restriction_indices]
  #z <- expand.grid(x,y)
  
  #gamma <- 2
  
  library(grid) 
  
  #df <- data.frame(x=diag(x)[list_complete$restriction_indices],y=diag(y)[list_complete$restriction_indices],dx=gamma*diag(list_complete$anisVector_1(x,y))[list_complete$restriction_indices],dy=gamma*diag(list_complete$anisVector_2(x,y))[list_complete$restriction_indices])
  df <- data.frame(x=loc.sub[,1],y=loc.sub[,2],dx=gamma*diag(list_complete$anisVector_1(x,y))[list_complete$restriction_indices],dy=gamma*diag(list_complete$anisVector_2(x,y))[list_complete$restriction_indices])
  #ggmap(map_white) + geom_segment(data=df,aes(x=x,y=y,xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.1,"cm")), color="red")
  return(ggmap(map_lite) + geom_segment(data=df,aes(x=x,y=y,xend=x+dx, yend=y+dy), arrow = arrow(length = unit(0.2,"cm")), color="black",size=size))
}

plot_vector_field(list_complete = list_complete)