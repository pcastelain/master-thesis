#####Completes the initial list with  the anisotropy functions
complete_list_anis <- function(initial_list, parameter, anis_mode = 1) {
  require(Matrix)
  
  if (length(parameter) < 2) {
    stop("Not enough parameters")
  }
  
  
  
  #Function to compute kappa and gamma (here called diagH)
  diagKappa <- function(x, y) {
    kappa <- as(x, "CsparseMatrix")
    
    
    diag(kappa) <- parameter[1]
    
    return(kappa)
  }
  
  diagH <- function(x, y) {
    gammaParam <- Matrix(0, dim(x)[1], dim(x)[2], sparse = TRUE)
    gammaParam[which((x != 0) & (y != 0))] <- parameter[2]
    
    return(gammaParam)
  }
  
  
  #Polynomial model
  if ((anis_mode == 1) && (length(parameter) > 2)) {
    anis_param <- parameter[-c(1:2)]
    
    
    
    #Here I should send warnings if things are not ok regarding parameter number
    
    #This vector format is v1 = a1 + b1 X + c1 Y + d1 XY + e1 X^2 + f1 y^2
    # v2 = a2 + b2 X + c2 Y + d2 XY + e2 X^2 + f2 y^2
    # and parameters are input in the following order (a1,a2,b1,b2,c1,c2,d1,d2,e1,e2,f1,f2)
    
    if (length(anis_param) < 2) {
      stop("Not enough Anisotropy parameters")
    } else if ((length(anis_param) %% 2) != 0) {
      stop("Not even number of Anisotropy parameters")
    }
    
    #We can model up to polynomial of degree 3 if there is less parameters we fill the end with 0
    if (length(anis_param) < 20) {
      anis_param[length(anis_param) + 1:20] <- 0
    }
    
    anisVector_1 <- function(x, y) {
      #With cosinus
      
      v1 <- Matrix(0, dim(x)[1], dim(x)[2], sparse = TRUE)
      
      indices_of_interest <- which((x != 0) & (y != 0))
      
      
      v1[indices_of_interest] <- anis_param[1] +
        anis_param[3] * x[indices_of_interest] +
        anis_param[5] * y[indices_of_interest] +
        anis_param[7] * x[indices_of_interest] * y[indices_of_interest] +
        anis_param[9] * x[indices_of_interest] ^ 2 +
        anis_param[11] * y[indices_of_interest] ^ 2 +
        anis_param[13] * x[indices_of_interest] ^ 2 * y[indices_of_interest] +
        anis_param[15] * y[indices_of_interest] ^ 2 * x[indices_of_interest] +
        anis_param[17] * x[indices_of_interest] ^ 3 +
        anis_param[19] * y[indices_of_interest] ^ 3
      
      return(v1)
    }
    
    anisVector_2 <- function(x, y) {
      # With sinus
      v2 <- Matrix(0, dim(x)[1], dim(x)[2], sparse = TRUE)
      
      indices_of_interest <- which((x != 0) & (y != 0))
      
      v2[indices_of_interest] <- anis_param[2] +
        anis_param[4] * x[indices_of_interest] +
        anis_param[6] * y[indices_of_interest] +
        anis_param[8] * x[indices_of_interest] * y[indices_of_interest] +
        anis_param[10] * x[indices_of_interest] ^ 2 +
        anis_param[12] * y[indices_of_interest] ^ 2 +
        anis_param[14] * x[indices_of_interest] ^ 2 * y[indices_of_interest] +
        anis_param[16] * y[indices_of_interest] ^ 2 * x[indices_of_interest] +
        anis_param[18] * x[indices_of_interest] ^ 3 +
        anis_param[20] * y[indices_of_interest] ^ 3
      return(v2)
    }
  } else if ((anis_mode == 2) && (length(parameter) > 2)) {
    ##################################################################
    ##                     Interpolation models                     ##
    ##################################################################
    
    require(akima)
    
    anis_param <- parameter[-c(1:2)]
    #This anisotropic function relates on bicubic splines models. We use a fixed grid and the parameters relates to
    
    #param order is X and then y !
    
    if (length(anis_param) < 4) {
      stop("Not enough Anisotropy parameters")
    } else if ((length(anis_param) %% 2) != 0) {
      stop("Not even number of Anisotropy parameters")
    }
    x_restricted <-
      diag(initial_list$D_i)[initial_list$restriction_indices]#,initial_list$restriction_indices])
    y_restricted <-
      diag(initial_list$D_j)[initial_list$restriction_indices]#,initial_list$restriction_indices])
    
    lenght_anis_grid <- sqrt(length(anis_param) / 2)
    #x_interp <- (1:lenght_anis_grid)/(lenght_anis_grid+1) * (max(x_restricted)-min(x_restricted)) + min(x_restricted)
    #y_interp <- (1:lenght_anis_grid)/(lenght_anis_grid+1) * (max(y_restricted)-min(y_restricted)) + min(y_restricted)
    x_interp <-
      seq(
        from = min(x_restricted),
        to = max(x_restricted),
        length.out = lenght_anis_grid
      )
    y_interp <-
      seq(
        from = min(y_restricted),
        to = max(y_restricted),
        length.out = lenght_anis_grid
      )
    
    
    anisVector_1 <- function(x, y) {
      v1 <- Matrix(0, dim(x)[1], dim(x)[2], sparse = TRUE)
      
      logic_vec <-
        (x != 0) * (y != 0)  #* (x < max(x_restricted)) * (x > min(x_restricted)) * (y < max(y_restricted)) * (y > min(y_restricted))
      indices_of_interest <-
        which(logic_vec == 1)#intersect(which(logic_vec==1),mat_close_indices) #Here it make anis function null outside certain location
      
      if (lenght_anis_grid < 4) {
        v1[indices_of_interest] <-
          akima::bilinear(
            x = x_interp,
            y = y_interp,
            z = matrix(anis_param[1:(length(anis_param) / 2)], ncol = lenght_anis_grid, nrow =
                         lenght_anis_grid) ,
            x0 = x[indices_of_interest],
            y0 = y[indices_of_interest]
          )$z
      } else{
        v1[indices_of_interest] <-
          akima::bicubic(
            x = x_interp,
            y = y_interp,
            z = matrix(anis_param[1:(length(anis_param) / 2)], ncol = lenght_anis_grid, nrow =
                         lenght_anis_grid) ,
            x0 = x[indices_of_interest],
            y0 = y[indices_of_interest]
          )$z
      }
      
      return(v1)
    }
    
    anisVector_2 <- function(x, y) {
      # With sinus
      v2 <- Matrix(0, dim(x)[1], dim(x)[2], sparse = TRUE)
      
      logic_vec <-
        (x != 0) * (y != 0) #* (x < max(x_restricted)) * (x > min(x_restricted)) * (y < max(y_restricted)) * (y > min(y_restricted))
      indices_of_interest <-
        which(logic_vec == 1) 
      
      if (lenght_anis_grid < 4) {
        v2[indices_of_interest] <-
          akima::bilinear(
            x = x_interp,
            y = y_interp,
            z = matrix(anis_param[(1 + length(anis_param) / 2):length(anis_param)], ncol =
                         lenght_anis_grid, nrow = lenght_anis_grid) ,
            x0 = x[indices_of_interest],
            y0 = y[indices_of_interest]
          )$z
        
      }
      else{
        v2[indices_of_interest] <-
          akima::bicubic(
            x = x_interp,
            y = y_interp,
            z = matrix(anis_param[(1 + length(anis_param) / 2):length(anis_param)], ncol =
                         lenght_anis_grid, nrow = lenght_anis_grid) ,
            x0 = x[indices_of_interest],
            y0 = y[indices_of_interest]
          )$z
      }
      
      
      return(v2)
    }
    
    
  }
  else{
    anisVector_1 <- function(x, y) {
      v1 <- Matrix(0, dim(x)[1], dim(x)[2], sparse = TRUE)
      return(v1)
    }
    
    anisVector_2 <- function(x, y) {
      v2 <- Matrix(0, dim(x)[1], dim(x)[2], sparse = TRUE)
      return(v2)
    }
  }
  
  ############## Il faut que je sorte le kappa de la liste pour le placer ici :)
  
  add_to_list <- list(anisVector_1, anisVector_2, diagH, diagKappa)
  names(add_to_list) <-
    c("anisVector_1", "anisVector_2", "diagH", "diagKappa")
  
  initial_list <- c(initial_list, add_to_list)
  
  return(initial_list)
}

##############Objective function that we optimize to find the anisotropy parameters
objective_function_gradient_score_anis = function(parameter,
                                                  excesses,
                                                  anis_list,
                                                  weighting_function,
                                                  weighting_function_derivative,
                                                  weight_coeff,
                                                  complete_list_anis = complete_list_anis,
                                                  nCores = 1,
                                                  r_func_vec,
                                                  cl = NULL,
                                                  compos_like_index_list = NULL,
                                                  first_param = NULL,
                                                  printPar = FALSE,
                                                  scoreFn = NULL,
                                                  anis_mode = NULL,
                                                  param_scale = 1) {
  # Print parameter to follow optimization evolution
  
  parameter <- parameter * param_scale
  
  if (printPar) {
    print(parameter)
  }
  
  if (!is.null(first_param)) {
    parameter <- c(first_param, parameter)
  }
  
  if (parameter[1] < 0 | parameter[2] < 0) {
    return(1e50)
  }
  
  if (is.null(anis_mode))
    anis_list <- complete_list_anis(anis_list, parameter = parameter)
  else{
    anis_list <-
      complete_list_anis(anis_list, parameter = parameter, anis_mode = anis_mode)
  }
  #Compute score
  if (is.null(scoreFn)) {
    scoreEstimation_v2(
      excesses,
      anis_list,
      weighting_function,
      weighting_function_derivative,
      weight_coeff = weight_coeff,
      r_func_vec = r_func_vec,
      nCores = nCores,
      cl = cl,
      compos_like_index_list = compos_like_index_list
    )
  }
  else{
    scoreFn(
      excesses,
      anis_list,
      weighting_function,
      weighting_function_derivative,
      weight_coeff = weight_coeff,
      r_func_vec = r_func_vec,
      nCores = nCores,
      cl = cl,
      compos_like_index_list = compos_like_index_list
    )
  }
  
}
