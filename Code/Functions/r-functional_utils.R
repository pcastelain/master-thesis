library(parallel)
library(foreach)
library(doParallel)

#Function to compute a linear-r functional provided the associated vector
linear_r <- function(database , r_functional_vec = NULL) {
  if (length(dim(database)) > 1){
    nb_sites_database <- length(database[1,])
  }
  else{
    nb_sites_database <- length(database)
  }

  
  if(is.null(r_functional_vec)){
    r_functional_vec <- rep(1/nb_sites_database,nb_sites_database)
  }
  
  if (nb_sites_database == length(r_functional_vec)){
    return(as.vector(database %*% r_functional_vec))
  }
  else{
    print("Error not same dimension")
  }
}

#Provide a vector based on a square region approach
generate_linear_r_vec <- function(coords){

  index <- as.logical(coords[,1] <= 5.25 & coords[,1] >= -1 & coords[,2] <= 52.5 & coords[,2] >= 48.75)

  # Here i should improve my integration process to take into account the earth ellipsoid ! It biased otherwise.
  return(index/sum(index))
}

transform_max_data <- function(data, time_window=12){
  
  registerDoParallel(detectCores())
  
  tot_time <- length(data[,1])
  
  max_data <- matrix(NA,nrow = tot_time, ncol = length(data[1,]))
  
  max_data <- foreach (i=c(1: tot_time), .combine=rbind) %dopar% {
    apply(data[max(1,i - time_window):min(tot_time,i+time_window),],MARGIN = 2 , FUN = max)
  }
  
  return(max_data)
  
  # 
  # for(i in c(1: tot_time)){
  #   if (i <= time_window){
  #     low <- 1
  #   }
  #   else{
  #     low <- i- time_window
  #   }
  #   if (i >= tot_time-time_window){
  #     high <- tot_time
  #   }
  #   else{
  #     high <- i + time_window
  #   }
  #   transform_max_data <- apply(data[low:high,],MARGIN = 2 , FUN = max)
  #   print(paste(i,"/",tot_time))
  # }
}

