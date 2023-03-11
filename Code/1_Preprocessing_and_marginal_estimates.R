##################################################################
##             Preprocessing and Marginal estimates             ##
##################################################################


# Load windgust dataset from ERA-Interim 
if (!file.exists("../Data/raw_database.RData")) download.file("https://renkulab.io/gitlab/paul.castelain/master-thesis/-/blob/master/Data/raw_database.RData","../Data/raw_database.RData")
load("../Data/raw_database.RData")
load("../Data/Land_sea_mask.RData")


#Load set of function related to r-functional
source("Functions/r-functional_utils.R")

# Load package for manipulation of coordinates
library(sp)
library(rgdal)
library(patchwork)
library(parallel)



#################Reordering indices 
reordering <- order(coords[,2],coords[,1])
coords <- coords[reordering,]
database <- database[,reordering]




reduce <- TRUE
#Reduce the size of the data `!`
if(reduce){
  
  upper <- c(10.5,57.75)#Area take from (de Fondeville et al., 2020)
  lower <- c(-8,44.25)
  
  index_spat <- (coords[,1] <= upper[1] & coords[,1] >= lower[1] & coords[,2] <= upper[2] & coords[,2] >= lower[2])
  
  coords <- coords[index_spat,]
  database <- database[,index_spat]
}

#A bolean mask, 1 if land, 0 if sea.
land_sea <- sapply(c(1:length(coords[,1])),function(i){Land_sea_mask[which(Land_sea_mask[,2] == coords[i,1] & Land_sea_mask[,3] == coords[i,2]),1]})
land_sea <- ifelse(land_sea>0.5,1,0)

##################################################################
##                          Footprints                          ##
##################################################################

time_window = 6 #number of 3hour timeframes before and after to use in footprint, 6*3=18h before and 18h after, 36hour footprint total

max_values_bol <- TRUE
max_data <- matrix(NA,nrow = nrow(database),ncol = ncol(database))
if(max_values_bol){
  max_data <- transform_max_data(database,time_window = time_window)
  rownames(max_data) <- NULL
}

# Create coordinates object
coords <- as.data.frame(coords)
names(coords) <-  c("X", "Y")

loc.sub <- cbind(coords$X,coords$Y)


#Load ggplot2 for plot
library(ggplot2)

# Set to true if you want to produce maps
map_plotting <- FALSE

if(map_plotting == TRUE){
  # Load package for map plotting
  library(ggmap)
  
  # Specify here your onwn google_api_key
  #api <- readLines("google_api_key")
  #register_google(key = api)
  
  # Downlaod map
  zoom_factor <- 1.1
  bbox <- c(min(coords$X)-(max(coords$X)-min(coords$X))*(zoom_factor-1),min(coords$Y)-(max(coords$Y)-min(coords$Y))*(zoom_factor-1),max(coords$X)+(max(coords$X)-min(coords$X))*(zoom_factor-1),max(coords$Y)+(max(coords$Y)-min(coords$Y))*(zoom_factor-1))
  
  #Here we use a boundary box to prevent the use of the zoom functionality (that needs google)
  map <- get_map(location = bbox, source="stamen", maptype = "toner" , color="bw")
  
  #Plot one sample from the dataset
  mapPlot(database[1,], loc.sub, title = "Raw wind data")
}

#################################################################
##             Risk definition and storm selection             ##
#################################################################


#Compute the r-functional linear vector
r_func_vec <- generate_linear_r_vec(loc.sub)
#un <- linear_r(bn,r_func_vec)

#n_r_func_vec 
n_r_func_vec <- r_func_vec*1/max(r_func_vec)#Normalized r_fun_vec for later use


#Is the risk based selection made on max_data or instantaneous! Instantaneous makes declustering more reliable.
instant <- TRUE
#Compute R(X) at each time sample
if(instant){
  r_x_t <- linear_r(database,r_func_vec)
} else{
  r_x_t <- linear_r(max_data,r_func_vec)
}

storm_index <- c() 

select_by_total_number <- 35 #(We are studying 35 years, 35 storms -> 1 big storm per year approx)
if (is.null(select_by_total_number)){
  find_max_r_x <- ifelse(r_x_t >= un,r_x_t,NA)
  
  tot_time <- length(find_max_r_x)
  while (sum(!is.na(find_max_r_x))!=0) {
    
    storm_index[length(storm_index)+1] <- which.max(find_max_r_x)
    
    find_max_r_x[max(storm_index[length(storm_index)]-2*time_window-1,1):min(tot_time,storm_index[length(storm_index)]+2*time_window+1)] <- rep(NA,1+min(tot_time - storm_index[length(storm_index)],2*time_window+1) + min(storm_index[length(storm_index)],1+2*time_window))
    
  }
}else{
  find_max_r_x <- r_x_t
  
  tot_time <- length(find_max_r_x)
  while (length(storm_index)<select_by_total_number) {
    
    storm_index[length(storm_index)+1] <- which.max(find_max_r_x)
    
    find_max_r_x[max(storm_index[length(storm_index)]-2*time_window-1,1):min(tot_time,storm_index[length(storm_index)]+2*time_window+1)] <- rep(NA,1+min(tot_time - storm_index[length(storm_index)],2*time_window+1) + min(storm_index[length(storm_index)],1+2*time_window))
    
  }
}


#define the new sample of storm footprints
storm_footprint <- max_data[storm_index,]

storm_dates <- dates[storm_index]


#################################################################
##            Define bn(s) based on storm selection            ##
#################################################################

quantileLevel_storm <- 0.5

storm_risk <- linear_r(storm_footprint,r_func_vec)

quantile_storm <- rep(NA, length(storm_footprint[1,]))


keep_looping <- TRUE
while(keep_looping){
  
  for(i in 1:length(storm_footprint[1,])){
    quantile_storm[i] <- quantile(storm_footprint[,i], quantileLevel_storm, na.rm = TRUE)
  }
  
  if (all(linear_r(quantile_storm,r_func_vec)<=storm_risk)){
    keep_looping = FALSE
  }else{
    quantileLevel_storm = quantileLevel_storm - 0.01
  }
}

bn <- quantile_storm

# un = r(bn)
un <- linear_r(bn,r_func_vec)
#Given by the explanation followin eq (6) in defondeville et al. 2020. We set r(A)=1 and r(B)=0 so bn_prime = un ! 
B <-  bn - un

##################################################################
##             Fit of pointwise GPD (not mandatory)             ##
##################################################################

#  for generalized Pareto distributions
scales <- rep(NA, length(bn))
shapes <- rep(NA, length(bn))

#List of fitted objects to check results
list_fit <- list()

# Local GP fit above bn(s) quantile
for(i in 1:length(scales)){
  print(paste(i, "/", length(scales)))
  
  # Likelihood based GPD fit
  fitGP <- evd::fpot(
    storm_footprint[, i],
    threshold = bn[i],
    std.err = FALSE,
    method = "Nelder-Mead",
    control = list(maxit = 10000)
  )
  
  # Store estimation results
  scales[i] <- fitGP$estimate[1]
  shapes[i] <- fitGP$estimate[2]
  list_fit[[i]] <- fitGP
}




#################################################################
##                 Fitting marginal parameters                 ##
#################################################################

#OUR loglike function does not depend on A as we decided A could be non continuous thus we define an and bn pointwise
#input X as MAtrix columns are loc rows are observations
indep_log <- function(x, xi, an, bn){
  
  one_event_log_like <- function(xj){

      if(sum(c(xj >= bn))>0){
        value <- (c(xj >= bn) %*% log(1/an * pmax(1 + xi * (xj-bn)/an,0)^(-(1/xi)-1)))# here removed the weight  /sum(c(xj >= bn))   ###########Here we need to check I dont believe eq. 23. It seems wrong compared to 21
      }
      else{
        value <- 0
      }
    
    # else{
    #   value <- (c(xj >= bn) %*% log(1/an * exp(-(xj-bn)/an))) #here removed the weight /sum(c(xj >= bn))
    # }
    #print(value)
    return(value)
  }
  if(is.vector(x)){
    vec_event_log <- one_event_log_like(x)
  }
  else{
    vec_event_log <- apply(x,1,one_event_log_like)
  }
 
  sum(vec_event_log)
}

#################################################################
##                  Grid search for xi and an                  ##
#################################################################
xi <- c()
an <- c()

xi_factor <- rep(0,length(loc.sub[,1]))

#xi_factor <- n_r_func_vec ##Set 2 xis one for the square of interest and one for the rest

for(xi_index in unique(xi_factor)){
  nCores <- detectCores()
  
  # Load package for parallelisation and create cluster instance
  if(nCores > 1){
    library(parallel)
    cl <- makeCluster(nCores,type="PSOCK") #Changed Fork to "PSOCK" as I'm running on windows
    clusterSetRNGStream(cl)
  } else {
    cl <- NULL
  }
  
  ############We split the locations for each cluster
  #index_loc <- c(1:length(bn))
  index_loc <- which(xi_factor == xi_index)
  list_index_per_cluster <- split(index_loc,((0:(length(index_loc)-1))%%nCores)+1)
  
  #grid_xi = seq(from= -0.05,to=-1,length.out = 200)
  
  grid_xi = seq(from= 1,to=-1,length.out = 200)
  
  cluster_fx_grid_search_xi_an <- function(vec_index){
    best_an <- matrix(NA, nrow = length(grid_xi),ncol = length(vec_index))
  
    for( k in c(1:length(grid_xi))){
      an_vec <- c()
      for(i in c(1:length(vec_index))){
        an_vec[i] <- suppressWarnings(tryCatch({
                              fit <- evd::fpot(storm_footprint[,vec_index[i]],bn[vec_index[i]],shape = grid_xi[k],lower=1,upper=30,method="Brent",std.err = FALSE)
                              fit$estimate[1]
                              }, error=function(e) { print(e); return(NULL) }))
        fit <- NULL
      }
      best_an[k,] <- an_vec
    }
    return(best_an)
  }
  
  xi_an_list_return_from_cluster <- NULL
  ###############Sending instruction to clusters
  clusterExport(cl,varlist=c("grid_xi","bn","storm_footprint"))  
  xi_an_list_return_from_cluster <- parLapply(cl,list_index_per_cluster,cluster_fx_grid_search_xi_an)
  
  stopCluster(cl)
  
  an_depending_on_xi <- matrix(NA,nrow = length(grid_xi),ncol = length(bn))
  for(i in c(1:nCores)){
    an_depending_on_xi[,list_index_per_cluster[[i]]] <- xi_an_list_return_from_cluster[[i]]
  }
  
  boxplot(t(an_depending_on_xi),names=grid_xi)
  
  value_for_xi <- c()
  for(i in c(1:length(grid_xi))){
    value_for_xi[i] <- -indep_log(x=storm_footprint[,index_loc],xi=grid_xi[i],an=an_depending_on_xi[i,index_loc],bn=bn[index_loc])
  }
  
  plot(grid_xi,value_for_xi)
  
  grid_xi[which.min(value_for_xi)]
  
  xi[index_loc] <- grid_xi[which.min(value_for_xi)]
  an[index_loc] <- an_depending_on_xi[which.min(value_for_xi),index_loc]
}


#################################################################
##                       Rescaling storm                       ##
#################################################################

rescaled_storm_footprint <- matrix(0,ncol=ncol(storm_footprint),nrow=nrow(storm_footprint))
for(i in c(1:nrow(storm_footprint))){
  rescaled_storm_footprint[i,] = pmax((1+xi*(storm_footprint[i,]-bn)/an),0)^(1/xi)
  #rescaled_storm_footprint[i,] = (1+xi*(storm_footprint[i,]-bn)/an)^(1/xi)
}

# Remove heavy and useless objects before saving workspace
rm("max_data","database")

# Saving workspace
save.image(file="../Tmp/Storm_Marginals.RData")


##################################################################
##                         Risk qq-plot                         ##
##################################################################

nb_sample <- 10000
confidence_inter <- 0.01

an_prime <- linear_r(an,r_func_vec)
bn_prime <- linear_r(bn,r_func_vec)

storm_risk_qq <- sort(linear_r(storm_footprint,r_func_vec))
nb_points <- length(storm_risk_qq)

#quantiles <- seq(from = 0, to= 1, length.out = nb_points+1)
#quantiles <- quantiles[-1] - quantiles[2]/2 #(1:n)/n+1

quantiles <- (1:nb_points)/(nb_points+1)

distrib_theory <- evd::qgpd(quantiles, loc=bn_prime, scale=an_prime, shape=unique(xi[n_r_func_vec==1]))

sample_mat <- t(sapply(c(1:nb_sample),function(k){sort(evd::rgpd(nb_points, loc=bn_prime, scale=an_prime, shape=unique(xi[n_r_func_vec==1])))}))

confidence_inter_vector <- sapply(c(1:nb_points), function(k){quantile(sample_mat[,k],probs=c(confidence_inter/2,1-(confidence_inter/2)))})
confidence_inter_low <- apply(sample_mat,2, function(k){quantile(k,probs=confidence_inter/2)})
confidence_inter_high <- apply(sample_mat,2, function(k){quantile(k,probs=1-(confidence_inter/2))})


lims <- c()
lims[1] <- min(storm_risk_qq,distrib_theory,confidence_inter_vector[1,],confidence_inter_vector[2,])
lims[2] <- max(storm_risk_qq,distrib_theory,confidence_inter_vector[1,],confidence_inter_vector[2,])


plot(distrib_theory,storm_risk_qq,xlim=lims,ylim=lims,pch=19,col="firebrick2",xlab="Model",ylab="Empirical")
abline(a=0,b=1)
lines(distrib_theory,confidence_inter_vector[1,],lty="dashed",lwd=2,col="skyblue4")
lines(distrib_theory,confidence_inter_vector[2,],lty="dashed",lwd=2,col="skyblue4")



##################################################################
##                       Marginal qq-plot                       ##
##################################################################


sample_index_big <- sample_index(loc.sub,5,4)

sample_index_small <- sample_index_big[c(1,5,8,18-5,21-5,25-5)]

layout(matrix(c(4,5,6,1,2,3), 2, 3, byrow = TRUE))

nb_sample <- 10000
confidence_inter <- 0.01

for(j in c(1:length(sample_index_small))){
  i <- sample_index_small[j]
  ordered_points <- sort(storm_footprint[,i])
  ordered_points <- ordered_points[which(ordered_points >= bn[i])] 
  nb_points <- length(ordered_points)
  
  quantiles <- (1:nb_points)/(nb_points+1)    
  
  
  distrib_theory <- evd::qgpd(quantiles, loc=bn[i], scale=an[i], shape=xi[i])
  
  
  sample_mat <- t(sapply(c(1:nb_sample),function(k){sort(evd::rgpd(nb_points, loc=bn[i], scale=an[i], shape=xi[i]))}))
  
  confidence_inter_vector <- sapply(c(1:nb_points), function(k){quantile(sample_mat[,k],probs=c(confidence_inter/2,1-(confidence_inter/2)))})
  
  quantile_points <- sapply(c(1:nb_points), function(k){ecdf(sample_mat[,k])(ordered_points[k])})
  
  qq_mat <- matrix(c(ordered_points,quantile_points,distrib_theory,confidence_inter_vector[1,],confidence_inter_vector[2,]),5,nb_points,byrow = TRUE)
  
  rownames(qq_mat) <- c("ordered_points","quantile_points","distrib_theory",paste(100*confidence_inter/2,"%",sep=""),paste(100*(1-confidence_inter/2),"%",sep=""))

  lims <- c()
  lims[1] <- min(qq_mat[3,],qq_mat[1,],qq_mat[4,],qq_mat[5,])
  lims[2] <- max(qq_mat[3,],qq_mat[1,],qq_mat[4,],qq_mat[5,])

  plot(qq_mat[3,],qq_mat[1,],xlim=lims,ylim=lims,pch=19,col="firebrick2",xlab="Model",ylab="Empirical", main=paste("Location",i))
  abline(a=0,b=1)
  lines(qq_mat[3,],qq_mat[4,],lty="dashed",lwd=2,col="skyblue4")
  lines(qq_mat[3,],qq_mat[5,],lty="dashed",lwd=2,col="skyblue4")

}
