# Function to estimate the empirical extremogram between two locations
computeExtremalCoeff <- function(pos1, pos2, database, local_quantiles) {
  theta <- length(which((database[, pos1] > local_quantiles[pos1] & !is.na(database[,pos1])) &
                          (database[, pos2] > local_quantiles[pos2] & !is.na(database[,pos2]))
  )) /
    length(which((database[, pos1] > local_quantiles[pos1] & !is.na(database[, pos1]))))
  return(ifelse(is.na(theta),0,theta))
}

#Create regularly spaced sample points. Returns indices !
sample_index <- function(loc,n_lon,n_lat){
  un_lon <- unique(loc[,1])
  un_lat <- unique(loc[,2])
  
  selected_lon <- (length(un_lon)+1)*(1:n_lon)/(n_lon+1)
  selected_lat <- (length(un_lat)+1)*(1:n_lat)/(n_lat+1)
  
  
  selected_lon <- ifelse(selected_lon > (length(un_lon)+1)/2,ceiling(selected_lon),floor(selected_lon))
  selected_lat <- ifelse(selected_lat > (length(un_lat)+1)/2,ceiling(selected_lat),floor(selected_lat))
  
  return(sort(which((loc[,1] %in% un_lon[selected_lon]) & (loc[,2] %in% un_lat[selected_lat]))))
}

#Return the covariance matrix based on a complete anisotropy list. The restricted parameter indicated wether to return the covariance for the whole SPDE domain
# or to have the covariance restricted to the studied locations in the middle of the domain
get_cov_mat <- function(complete_list,restricted=TRUE){
  
  attach(complete_list)
  H_11 <- compute_A_H_11(coordinates_H_11_x_right,
                         coordinates_H_11_x_left,
                         coordinates_H_11_y,
                         sign_H_11,
                         anisVector_1,
                         diagH)
  H_22 <- compute_A_H_22(coordinates_H_22_x,
                         coordinates_H_22_y_top,
                         coordinates_H_22_y_bot,
                         sign_H_22,
                         anisVector_2,
                         diagH)
  H_12 <- compute_A_H_12(coordinates_H_12_x,
                         coordinates_H_12_y_top,
                         coordinates_H_12_y_bot,
                         sign_H_12_top,
                         sign_H_12_bot,
                         anisVector_1,
                         anisVector_2,
                         diagH)
  H_21 <- compute_A_H_21(coordinates_H_21_x_right,
                         coordinates_H_21_x_left,
                         coordinates_H_21_y,
                         sign_H_21_right,
                         sign_H_21_left,
                         anisVector_1,
                         anisVector_2,
                         diagH)
  
  ######## COMPUTE PRECISION MATRIX #########
  A_H <- (h_y / h_x) * H_11 + (h_x / h_y) * H_22 + 1 / 4  * (H_12 + H_21)  # i don't get why there are h_x,h_y coefficients
  matrix_A <- D_V %*% diagKappa(D_i,D_j) - A_H #D_V = I * V with V=h_x*h_y, its the area of each square
  precisionMatrix_Q <- Matrix::t(matrix_A) %*% D_V_inv %*% matrix_A 
  
  covarianceMatrix  <-  solve(precisionMatrix_Q) 
  
  if(! is.null(restriction_indices) & restricted){
    covarianceMatrix <- covarianceMatrix[restriction_indices,restriction_indices]
  }
  
  detach(complete_list)
  return(covarianceMatrix)
}


generate_sample <- function(complete_list,xi,an,bn,return_level_years=NULL,loc_ref=1){
  require(EnvStats)
  
  covarianceMatrix <- get_cov_mat(complete_list)
  
  #loc_ref <- 1
  
  risk <- function(vec){
    linear_r(vec,r_func_vec)
  }
  
  xi_risk <- sum(xi* r_func_vec) #The xi is supposed to be stationary over the risk region
  
  A <- an/risk(an)
  
  ##################I need to define u>0 st it does maximise the success of the R(Ay^xi)>1 while including the whole subset
  ## I can inverse the functions
  risk_A <- r_func_vec * A
  
  risk_A <- risk_A[which(risk_A!=0)]
  
  xi_risk_A <- xi[which(risk_A!=0)]
  
  minimum_norm_1 <- sum((1/risk_A)^(1/(xi_risk_A-1)))^((xi_risk-1)/xi_risk)
  
  u <- minimum_norm_1         
  
  Yr <- rep(0,length(an))
  
  ##########We added a condition on Yr ==0 cause it was exploding the risk otherwise with a division by 0 when the xi param is negative
  while(any(Yr == 0) || (risk(A*(Yr)^xi) > 1)){##################WARNING Repasser pour nettoyer !
    
    #R1 is unit pareto, density (1-1/r)
    R1 <- rpareto(1,1,1) 
    
    #wE GENERATE W1 log gaussian according to the covariance found 
    sim <- MASS::mvrnorm(n = 1, rep(0, length(an)), covarianceMatrix) 
    W1 <- sim - sim[loc_ref]
    W1[-loc_ref] <- W1[-loc_ref] - diag(covarianceMatrix)[-loc_ref] / 2
    W1 <- exp(W1)/norm(as.matrix(exp(W1)),type="1")
    
    Yr = R1*W1*u
    print(paste("Risk is: ",risk(A*(Yr)^xi)," | Norm1 of Yr :",norm(as.matrix(Yr),type="1"), " | Norm1 of W1: ", norm(as.matrix(W1),type="1")))
  }
  W2 <- (A*(Yr)^xi)/norm(as.matrix((A*(Yr)^xi)),type="2")#####En fait ca va#######J'ai un soucis avec la norme ! Il faut que je pondère en fonction de projected 
  
  #####If i fix R2 to a return level i'm interested in I'll get samples with the good return level
  if(is.null(return_level_years)){
    R2 <- rpareto(1,1,1) 
  }
  else{
    R2 <- qpareto(p = (1-1/return_level_years),location = 1, shape = 1)
  }
  
  
  P <- risk(an)*(xi^(-1))*R2^(xi)*W2/risk(W2)+ bn - (xi^(-1))*an
  return(P)
  #nice_plot(P,loc = loc.sub,map=map)
}

nice_plot <- function(dataVec,loc,map=NULL,alpha=NULL,lims = NULL,main=NULL,pal=rev(RColorBrewer::brewer.pal(7,"RdYlBu")),legend.title=NULL){
  require(ggmap)
  
  if(is.null(lims)) lims <- c(min(dataVec),max(dataVec))
  
  h_x <- abs(unique(loc[,1])[2]-unique(loc[,1])[1])
  h_y <- abs(unique(loc[,2])[2]-unique(loc[,2])[1])
  
  dataToPlot = data.frame(X = rep(1, 4 * length(dataVec)), Y = rep(1, 4 * length(dataVec)), Value = rep(1, 4 * length(dataVec)), Group = rep(1, 4 * length(dataVec)))
  for(i in 1:length(dataVec)){
    dataToPlot[(i - 1) * 4 + 1,] <- c(loc[i,1] + h_x / 2, loc[i,2] + h_y / 2, dataVec[i], i)
    dataToPlot[(i - 1) * 4 + 2,] <- c(loc[i,1] + h_x / 2 , loc[i,2] - h_y / 2, dataVec[i], i)
    dataToPlot[(i - 1) * 4 + 3,] <- c(loc[i,1] - h_x / 2, loc[i,2] - h_y / 2, dataVec[i] ,i)
    dataToPlot[(i - 1) * 4 + 4,] <- c(loc[i,1] - h_x / 2, loc[i,2] + h_y / 2, dataVec[i], i)
    
  }
  
  if(is.null(map)){
    if(is.null(alpha)){
      alpha <- 0.9
    }
    mapPoints <- ggplot(dataToPlot, aes(x = x, y = y)) + geom_polygon(data = dataToPlot, aes(x=X, y = Y,group = Group, fill = Value), alpha = alpha, color = NA) +
      scale_fill_gradientn(colors = pal,limits=lims, name = legend.title) +  theme(plot.title = element_text(hjust = 0.5))
  }
  else{
    if(is.null(alpha)){
      alpha <- 0.95 #0.7 in stamen terrrain bw, for stamen toner 0.95 is better
    }
    mapPoints <-  ggmap(map) + geom_polygon(data = dataToPlot, aes(x=X, y = Y,group = Group, fill = Value), alpha = alpha, color = NA) +
      scale_fill_gradientn(colors = pal,limits=lims, name= legend.title) +  theme(plot.title = element_text(hjust = 0.5))
  }
  if(is.null(legend.title)){
    mapPoints <- mapPoints + theme(legend.title = element_blank())
  }
  
  if(!is.null(main)){
    mapPoints <- mapPoints + labs(title = main, fill = "")
  }
  mapPoints
}

nice_layer <- function(loc,color_fill,alpha_fill=0.95,standard_h = TRUE){
  
  if (is.vector(loc)){
    length_vec <- 1
    loc <-  rbind(loc,loc)
  }
  else{
    length_vec <- length(loc[,1])
  }
  
  
  if(standard_h == FALSE){
    
    if(length_vec == 1){
      loc <- t(as.matrix(loc))
      h_x <- 0.75#This is in case there is only one point ! then the h_x=0=h_y
      h_y <- 0.75
    }
    else{
      h_x <- abs(unique(loc[,1])[2]-unique(loc[,1])[1])
      h_y <- abs(unique(loc[,2])[2]-unique(loc[,2])[1])
    }
    
    if(h_x == 0) h_x <- 0.75#This is in case there is only one point ! then the h_x=0=h_y
    if(h_y == 0) h_y <- 0.75
    
  }else{
    h_x <- 0.75
    h_y <- 0.75
  }
  
  dataToPlot = data.frame(X = rep(1, 4 * length_vec), Y = rep(1, 4 * length_vec), Group = rep(1, 4 * length_vec))
  for(i in 1:length_vec){
    dataToPlot[(i - 1) * 4 + 1,] <- c(loc[i,1] + h_x / 2, loc[i,2] + h_y / 2, i)
    dataToPlot[(i - 1) * 4 + 2,] <- c(loc[i,1] + h_x / 2 , loc[i,2] - h_y / 2, i)
    dataToPlot[(i - 1) * 4 + 3,] <- c(loc[i,1] - h_x / 2, loc[i,2] - h_y / 2 ,i)
    dataToPlot[(i - 1) * 4 + 4,] <- c(loc[i,1] - h_x / 2, loc[i,2] + h_y / 2, i)
    
  }
  return(geom_polygon(data = dataToPlot, aes(x=X, y = Y,group = Group), fill = color_fill, alpha = alpha_fill, color = NA))
}
varying_layer <- function(dataVec,loc,alpha,map_to_plot=map){
  
  h_x <- abs(unique(loc[,1])[2]-unique(loc[,1])[1])
  h_y <- abs(unique(loc[,2])[2]-unique(loc[,2])[1])
  
  dataToPlot = data.frame(X = rep(1, 4 * length(dataVec)), Y = rep(1, 4 * length(dataVec)), Value = rep(1, 4 * length(dataVec)), Group = rep(1, 4 * length(dataVec)))
  for(i in 1:length(dataVec)){
    dataToPlot[(i - 1) * 4 + 1,] <- c(loc[i,1] + h_x / 2, loc[i,2] + h_y / 2, dataVec[i], i)
    dataToPlot[(i - 1) * 4 + 2,] <- c(loc[i,1] + h_x / 2 , loc[i,2] - h_y / 2, dataVec[i], i)
    dataToPlot[(i - 1) * 4 + 3,] <- c(loc[i,1] - h_x / 2, loc[i,2] - h_y / 2, dataVec[i] ,i)
    dataToPlot[(i - 1) * 4 + 4,] <- c(loc[i,1] - h_x / 2, loc[i,2] + h_y / 2, dataVec[i], i)
    
  }
  
  mapPoints <- ggmap(map_to_plot) + geom_polygon(data = dataToPlot, aes(x=X, y = Y,group = Group, fill = Value), alpha = alpha, color = NA) +
    scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(7,"RdYlBu")),limits = c(0,1))
  return(mapPoints)
}