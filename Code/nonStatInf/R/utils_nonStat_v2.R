###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ######
###### Functions for non-stationary Matern field ######
###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ######

### INITialise, return dataframe of coordinates ###

#As input: Loc that has the format expand grid (M,N)
#Usage:
# init_anis(M=length(unique(loc.sub[,1])),N=length(unique(loc.sub[,2])),loc=as.data.frame(projected_coordinates))

#' @export
init_anis <- function(M,N,loc,restriction_indices=NULL,land_sea_map=NULL){
  #Library
  library(Matrix)
  library(ggplot2)

  if (dim(loc)[1] != (M*N)){
    stop("`loc` dimensions not matching with `N` and `M`.")
  }

  names(loc) <- c("X","Y")

  h_x <- c()
  h_y <- c()

  #Define h_x and h_y for each posRef
  for(i in 0:(M -1)){
    for(j in (0:(N - 1))){

      posRef <- j * M + i + 1

      if(i==0){
        i_p <- i #previous
        i_n <- i + 1 #next
        alpha_i = 1
      }
      else if(i==(M-1)){
        i_p <- i - 1 #previous
        i_n <- i #next
        alpha_i = 1
      }
      else{
        i_p <- i - 1 #previous
        i_n <- i + 1 #next
        alpha_i = 2
      }

      if(j==0){
        j_p <- j #previous
        j_n <- j + 1 #next
        alpha_j = 1
      }
      else if(j==(N-1)){
        j_p <- j - 1 #previous
        j_n <- j #next
        alpha_j = 1
      }
      else{
        j_p <- j - 1 #previous
        j_n <- j + 1 #next
        alpha_j = 2
      }

      h_x[posRef] <- (loc$X[j * M + i_n + 1]-loc$X[j * M + i_p + 1])/alpha_i
      h_y[posRef] <- (loc$Y[j_n * M + i + 1]-loc$Y[j_p * M + i + 1])/alpha_j

    }
  }

  V <- h_x * h_y

  #Anisotropy parameters


  #Define areal matrix
  D_V <- matrix(0, M*N, M*N)
  D_V_inv <- matrix(0, M*N, M*N)
  diag(D_V) <- V
  diag(D_V_inv) <- 1/V
  D_V <- as(D_V, "CsparseMatrix") #Change matrix to sparse matrix
  D_V_inv <- as(D_V_inv, "CsparseMatrix")

  #Diagonal coordinate matrix for calculation of varying kappa
  D_i <- matrix(0, M*N, M*N)
  diag(D_i) <- loc$X
  D_i <- as(D_i, "CsparseMatrix")

  D_j <- matrix(0, M*N, M*N)
  diag(D_j) <- loc$Y
  D_j <- as(D_j, "CsparseMatrix")

  #Function for calculation kappa
  kappa <- function(A,B){
    kappa <- A
    diag(kappa) <- 0.1
    kappa
  }

  #Diagonal matrix for kappa
  D_kappa <- kappa(D_i,D_j)

  #Matrix representation for A
  coordinates_H_11_x_right <- as(matrix(0, M * N, M * N), "CsparseMatrix")
  coordinates_H_11_x_left <- as(matrix(0, M * N, M * N), "CsparseMatrix")
  coordinates_H_11_y <- as(matrix(0, M * N, M * N), "CsparseMatrix")

  coordinates_H_22_x <- as(matrix(0, M * N, M * N), "CsparseMatrix")
  coordinates_H_22_y_top <- as(matrix(0, M * N, M * N), "CsparseMatrix")
  coordinates_H_22_y_bot <- as(matrix(0, M * N, M * N), "CsparseMatrix")

  coordinates_H_12_x <- as(matrix(0, M * N, M * N), "CsparseMatrix")
  coordinates_H_12_y_top <- as(matrix(0, M * N, M * N), "CsparseMatrix")
  coordinates_H_12_y_bot <- as(matrix(0, M * N, M * N), "CsparseMatrix")

  coordinates_H_21_x_right <- as(matrix(0, M * N, M * N), "CsparseMatrix")
  coordinates_H_21_x_left <- as(matrix(0, M * N, M * N), "CsparseMatrix")
  coordinates_H_21_y <- as(matrix(0, M * N, M * N), "CsparseMatrix")

  sign_H_11 <- as(matrix(0, M * N, M * N), "CsparseMatrix")
  sign_H_22 <- as(matrix(0, M * N, M * N), "CsparseMatrix")
  sign_H_21_left <- as(matrix(0, M * N, M * N), "CsparseMatrix")
  sign_H_21_right <- as(matrix(0, M * N, M * N), "CsparseMatrix")
  sign_H_12_top <- as(matrix(0, M * N, M * N), "CsparseMatrix")
  sign_H_12_bot <- as(matrix(0, M * N, M * N), "CsparseMatrix")

  for(i in 0:(M -1)){
    for(j in (0:(N - 1))){

      posRef <- j * M + i + 1

      coordinates_H_11_x_right[posRef,posRef] <- loc$X[posRef] + h_x[posRef] / 2
      coordinates_H_11_x_left[posRef,posRef] <- loc$X[posRef] - h_x[posRef] / 2
      coordinates_H_11_y[posRef,posRef] <- loc$Y[posRef]

      sign_H_11[posRef,posRef] <- -1

      coordinates_H_22_y_top[posRef,posRef] <- loc$Y[posRef] + h_y[posRef] / 2
      coordinates_H_22_y_bot[posRef,posRef] <- loc$Y[posRef] - h_y[posRef] / 2
      coordinates_H_22_x[posRef,posRef] <- loc$X[posRef]

      sign_H_22[posRef,posRef] <- -1

      i_p <- i - 1 #previous
      i_n <- i + 1 #next

      j_p <- j - 1
      j_n <- j + 1

      if(i_p > -1){
        #### USE RIGHT FOR BOTH IN ORDER TO SYMMETRIZE #### #pas compris
        coordinates_H_11_x_left[posRef, j * M + i_p + 1] <- loc$X[posRef] - h_x[posRef] / 2
        coordinates_H_11_y[posRef,j * M + i_p + 1] <- loc$Y[posRef]
        sign_H_11[posRef,j * M + i_p + 1] <- 1

        coordinates_H_12_y_top[posRef, j * M + i_p + 1] <- loc$Y[posRef] + h_y[posRef] / 2
        coordinates_H_12_y_bot[posRef, j * M + i_p + 1] <- loc$Y[posRef] - h_y[posRef] / 2
        coordinates_H_12_x[posRef,j * M + i_p + 1] <- loc$X[posRef]
        sign_H_12_top[posRef,j * M + i_p + 1] <- -1
        sign_H_12_bot[posRef,j * M + i_p + 1] <- 1
      }

      if(i_n < M){
        #### USE RIGHT FOR BOTH IN ORDER TO SYMMETRIZE ####
        coordinates_H_11_x_right[posRef, j * M + i_n + 1] <- loc$X[posRef] + h_x[posRef] / 2
        coordinates_H_11_y[posRef,j * M + i_n + 1] <- loc$Y[posRef]
        sign_H_11[posRef,j * M + i_n + 1] <- 1

        coordinates_H_12_y_top[posRef, j * M + i_n + 1] <- loc$Y[posRef] + h_y[posRef] / 2
        coordinates_H_12_y_bot[posRef, j * M + i_n + 1] <- loc$Y[posRef] - h_y[posRef] / 2
        coordinates_H_12_x[posRef,j * M + i_n + 1] <- loc$X[posRef]
        sign_H_12_top[posRef,j * M + i_n + 1] <- 1
        sign_H_12_bot[posRef,j * M + i_n + 1] <- -1
      }

      if(j_p > -1){
        #### USE TOP FOR BOTH IN ORDER TO SYMMETRIZE ####
        coordinates_H_22_y_top[posRef,j_p * M + i + 1] <- loc$Y[posRef] - h_y[posRef] / 2
        coordinates_H_22_x[posRef,j_p * M + i + 1] <- loc$X[posRef]
        sign_H_22[posRef,j_p * M + i + 1] <- 1

        coordinates_H_21_x_right[posRef, j_p * M + i + 1] <- loc$X[posRef] + h_x[posRef] / 2
        coordinates_H_21_x_left[posRef, j_p * M + i + 1] <- loc$X[posRef] - h_x[posRef] / 2
        coordinates_H_21_y[posRef, j_p * M + i + 1] <- loc$Y[posRef]
        sign_H_21_right[posRef, j_p * M + i + 1] <- 1
        sign_H_21_left[posRef, j_p * M + i + 1] <- -1
      }

      if(j_n < N){
        #### USE TOP FOR BOTH IN ORDER TO SYMMETRIZE ####
        coordinates_H_22_y_bot[posRef,j_n * M + i + 1] <- loc$Y[posRef] + h_y[posRef] / 2
        coordinates_H_22_x[posRef,j_n * M + i + 1] <- loc$X[posRef]
        sign_H_22[posRef,j_n * M + i + 1] <- 1

        coordinates_H_21_x_right[posRef, j_n * M + i + 1] <- loc$X[posRef] + h_x[posRef] / 2
        coordinates_H_21_x_left[posRef, j_n * M + i + 1] <- loc$X[posRef] - h_x[posRef] / 2
        coordinates_H_21_y[posRef, j_n * M + i + 1] <- loc$Y[posRef]
        sign_H_21_right[posRef, j_n * M + i + 1] <- -1
        sign_H_21_left[posRef, j_n * M + i + 1] <- 1
      }

      if((i_p > -1) & (j_p > -1)){
        coordinates_H_12_y_bot[posRef, j_p * M + i_p + 1] <- loc$Y[posRef] - h_y[posRef] / 2
        coordinates_H_12_x[posRef,j_p * M + i_p + 1] <- loc$X[posRef]
        sign_H_12_bot[posRef,j_p * M + i_p + 1] <- 1

        coordinates_H_21_x_left[posRef, j_p * M + i_p + 1] <- loc$X[posRef] - h_x[posRef] / 2
        coordinates_H_21_y[posRef, j_p * M + i_p + 1] <- loc$Y[posRef]
        sign_H_21_left[posRef, j_p * M + i_p + 1] <- 1
      }

      if((i_n < M) & (j_p > -1)){
        coordinates_H_12_y_bot[posRef, j_p * M + i_n + 1] <- loc$Y[posRef] - h_y[posRef] / 2
        coordinates_H_12_x[posRef, j_p * M + i_n + 1] <- loc$X[posRef]
        sign_H_12_bot[posRef,j_p * M + i_n + 1] <- -1

        coordinates_H_21_x_right[posRef, j_p * M + i_n + 1] <- loc$X[posRef] + h_x[posRef] / 2
        coordinates_H_21_y[posRef, j_p * M + i_n + 1] <- loc$Y[posRef]
        sign_H_21_right[posRef, j_p * M + i_n + 1] <- -1
      }

      if((i_p > -1) & (j_n < N)){
        coordinates_H_12_y_top[posRef, j_n * M + i_p + 1] <- loc$Y[posRef] + h_y[posRef] / 2
        coordinates_H_12_x[posRef,j_n * M + i_p + 1] <- loc$X[posRef]
        sign_H_12_top[posRef,j_n * M + i_p + 1] <- -1

        coordinates_H_21_x_left[posRef, j_n * M + i_p + 1] <- loc$X[posRef] - h_x[posRef] / 2
        coordinates_H_21_y[posRef, j_n * M + i_p + 1] <- loc$Y[posRef]
        sign_H_21_left[posRef, j_n * M + i_p + 1] <- -1
      }

      if((i_n < M) & (j_n < N)){
        coordinates_H_12_y_top[posRef, j_n * M + i_n + 1] <- loc$Y[posRef] + h_y[posRef] / 2
        coordinates_H_12_x[posRef,j_n * M + i_n + 1] <- loc$X[posRef]
        sign_H_12_top[posRef,j_n * M + i_n + 1] <- 1

        coordinates_H_21_x_right[posRef, j_n * M + i_n + 1] <- loc$X[posRef] + h_x[posRef] / 2
        coordinates_H_21_y[posRef, j_n * M + i_n + 1] <- loc$Y[posRef]
        sign_H_21_right[posRef, j_n * M + i_n + 1] <- 1
      }
    }
  }

  if(is.null(restriction_indices)){
      restriction_indices <- c(1:(M*N))
  }

  #Diagonal matrix for kappa
  return_list <- list(h_x,h_y,D_V, D_V_inv,D_i,D_j,
                      coordinates_H_11_x_right,coordinates_H_11_x_left,coordinates_H_11_y,
                      coordinates_H_22_x,coordinates_H_22_y_top,coordinates_H_22_y_bot,
                      coordinates_H_12_x,coordinates_H_12_y_top,coordinates_H_12_y_bot,
                      coordinates_H_21_x_right,coordinates_H_21_x_left,coordinates_H_21_y,
                      sign_H_11,sign_H_22,sign_H_21_left,sign_H_21_right,sign_H_12_top,sign_H_12_bot,
                      restriction_indices,land_sea_map)
  names(return_list) <- c("h_x","h_y","D_V","D_V_inv","D_i","D_j",
                          "coordinates_H_11_x_right","coordinates_H_11_x_left","coordinates_H_11_y",
                          "coordinates_H_22_x","coordinates_H_22_y_top","coordinates_H_22_y_bot",
                          "coordinates_H_12_x","coordinates_H_12_y_top","coordinates_H_12_y_bot",
                          "coordinates_H_21_x_right","coordinates_H_21_x_left","coordinates_H_21_y",
                          "sign_H_11","sign_H_22","sign_H_21_left","sign_H_21_right","sign_H_12_top","sign_H_12_bot",
                          "restriction_indices","land_sea_map")
  return(return_list)
}


# ### ANISOTROPIC VECTOR FIELD  DEFINITION ###
# #Maybe need to rewrite as we are not using full potential of matrix sparcity because we call X and Y indices using logical values of (x>0,y>0)
# anisVector_1 <- function(x,y){ #With cosinus
#   v1 <- Matrix(0, dim(x)[1], dim(x)[2], sparse = TRUE)
#
#   dist_refPoint <- sqrt((x[(x > 0) & (y > 0)] - 20)^2 +  (y[(x > 0) & (y > 0)] - 20) ^2)
#   #angles <- atan((y[(x > 0) & (y > 0)] - 20) / (x[(x > 0) & (y > 0)] - 20))
#
#   v1[(x > 0) & (y > 0)] <- pmax(0,(40 - dist_refPoint)) * cos(dist_refPoint / 5)
#   # v1[(x > 0) & (y > 0)] <- 1.5
#   return(v1)
# }
#
# anisVector_2 <- function(x,y){ # With sinus
#   v2 <- Matrix(0, dim(x)[1], dim(x)[2], sparse = TRUE)
#
#   dist_refPoint <- sqrt((x[(x > 0) & (y > 0)] - 20)^2 +  (y[(x > 0) & (y > 0)] - 20) ^2)
#
#   v2[(x > 0) & (y > 0)] <-  pmax(0,(40 - dist_refPoint)) * sin(dist_refPoint / 5)
#   # v2[(x > 0) & (y > 0)] <-  1.5
#   return(v2)
# }
#
# diagH <- function(x,y){
#   gammaParam <- Matrix(0, dim(x)[1], dim(x)[2], sparse = TRUE)
#   gammaParam[(x > 0) & (y > 0)] <- 1 #1.2
#   return(gammaParam)
# }


#' @export
compute_A_H_11 <- function(coordinates_H_11_x_right,
                           coordinates_H_11_x_left,
                           coordinates_H_11_y,
                           sign_H_11,
                           anisVector_1,
                           diagH){
  anisVec_right_1 <- anisVector_1(coordinates_H_11_x_right, coordinates_H_11_y)
  anisVec_left_1 <- anisVector_1(coordinates_H_11_x_left, coordinates_H_11_y)
  gammaParam_right <- diagH(coordinates_H_11_x_right, coordinates_H_11_y)
  gammaParam_left <- diagH(coordinates_H_11_x_left, coordinates_H_11_y)

  H_11_left <-  gammaParam_left + anisVec_left_1^2 #It is squared because it is v*vT
  H_11_right <-  gammaParam_right + anisVec_right_1^2

  H_11 <- sign_H_11 * (H_11_left + H_11_right)

  return(H_11)
}
#' @export
compute_A_H_22 <- function(coordinates_H_22_x,
                           coordinates_H_22_y_top,
                           coordinates_H_22_y_bot,
                           sign_H_22,
                           anisVector_2,
                           diagH){
  anisVec_top_2 <- anisVector_2(coordinates_H_22_x, coordinates_H_22_y_top)
  anisVec_bot_2 <- anisVector_2(coordinates_H_22_x, coordinates_H_22_y_bot)
  gammaParam_top <- diagH(coordinates_H_22_x, coordinates_H_22_y_top)
  gammaParam_bot <- diagH(coordinates_H_22_x, coordinates_H_22_y_bot)

  H_22 <-  sign_H_22 * ((gammaParam_top + anisVec_top_2^2) + (gammaParam_bot +  anisVec_bot_2^2))

  return(H_22)
}
#' @export
compute_A_H_12 <- function(coordinates_H_12_x,
                           coordinates_H_12_y_top,
                           coordinates_H_12_y_bot,
                           sign_H_12_top,
                           sign_H_12_bot,
                           anisVector_1,
                           anisVector_2,
                           diagH){
  anisVec_top_1 <- anisVector_1(coordinates_H_12_x, coordinates_H_12_y_top)
  anisVec_top_2 <- anisVector_2(coordinates_H_12_x, coordinates_H_12_y_top)
  anisVec_bot_1 <- anisVector_1(coordinates_H_12_x, coordinates_H_12_y_bot)
  anisVec_bot_2 <- anisVector_2(coordinates_H_12_x, coordinates_H_12_y_bot)
  gammaParam_top <- diagH(coordinates_H_12_x, coordinates_H_12_y_top)
  gammaParam_bot <- diagH(coordinates_H_12_x, coordinates_H_12_y_bot)

  H_12 <- sign_H_12_top *  (anisVec_top_2 * anisVec_top_1) +
    sign_H_12_bot * (anisVec_bot_2 * anisVec_bot_1)

  return(H_12)
}
#' @export
compute_A_H_21 <- function(coordinates_H_21_x_right,
                           coordinates_H_21_x_left,
                           coordinates_H_21_y,
                           sign_H_21_right,
                           sign_H_21_left,
                           anisVector_1,
                           anisVector_2,
                           diagH){
  anisVec_right_1 <- anisVector_1(coordinates_H_21_x_right, coordinates_H_21_y)
  anisVec_right_2 <- anisVector_2(coordinates_H_21_x_right, coordinates_H_21_y)
  anisVec_left_1 <- anisVector_1(coordinates_H_21_x_left, coordinates_H_21_y)
  anisVec_left_2 <- anisVector_2(coordinates_H_21_x_left, coordinates_H_21_y)
  gammaParam_right <- diagH(coordinates_H_21_x_right, coordinates_H_21_y)
  gammaParam_Left <- diagH(coordinates_H_21_x_left, coordinates_H_21_y)

  H_21 <-  sign_H_21_right * (anisVec_right_1 * anisVec_right_2) +
    sign_H_21_left * (anisVec_left_1 * anisVec_left_2)

  return(H_21)
}



