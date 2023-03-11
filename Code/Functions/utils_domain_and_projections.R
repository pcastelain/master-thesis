###Returns projected coordinates
loc_to_proj <- function(loc,scaling=FALSE,scale=NULL,center=NULL,projecting=FALSE,fixed_lat=NULL){
  ############Nudge effect !
  #If we have x or y values equal to 0, the code will break ! To prevent this we slightly move the centering values !
  nudge <- 0.01
  
  require(sp)
  require(rgdal)
  # Create coordinates object
  loc <- as.data.frame(loc)
  names(loc) <-  c("X", "Y")
  coordinates(loc) <- c("X", "Y")
  
  # Define system of coordinates
  proj4string(loc) <- CRS("+proj=longlat +datum=WGS84")
  
  
  if(!is.null(fixed_lat)){
    values_to_return <- 6371 * cbind(2*pi*loc$X/360*cos(2*pi*fixed_lat/360),2*pi*loc$Y/360)
    values_to_return <- scale(values_to_return,center=TRUE,scale=FALSE)
    if(any(abs(values_to_return[,1])<nudge | abs(values_to_return[,1]) < nudge)){
      attr(values_to_return,"scaled:center") <- attr(values_to_return,"scaled:center") + c(-nudge,-nudge)
      values_to_return <- values_to_return + nudge
      warning("Some x or y were to equal 0 ! Scaling has been adjusted to correct this")
    }
    return(values_to_return)
  }
  else{
    
    if(projecting){
      # Project lon/lat coordinates to UTM (m) system of coordinates
      loc_proj <-  spTransform(loc, CRS("+proj=utm +zone=31 +datum=WGS84")) #The zone is based on location
      
      # Rescale to have new coordinates in km
      loc_proj@coords <-  loc_proj@coords / 1000
    }else{
      loc_proj <- loc
    }
    
    
    
    if(scaling){
      if(!is.null(scale) | !is.null(center)){
        if(is.null(scale)){
          scale <- c(1,1)
        }
        if(is.null(center)){
          center <- c(0,0)
        }
        loc_proj@coords <- scale(loc_proj@coords,center = center,scale=scale)
        if(any(loc_proj@coords[,1]== 0 | loc_proj@coords[,2] == 0)){
          warning("Some x or y values equal 0 ! This will make inference fail later !")
        }
      }
      else{
        loc_proj@coords <- scale(loc_proj@coords)
        if(any(loc_proj@coords[,1]== 0 | loc_proj@coords[,2] == 0)){
          attr(loc_proj@coords,"scaled:center") <- attr(loc_proj@coords,"scaled:center") + c(-nudge,-nudge)
          loc_proj@coords <- loc_proj@coords + nudge
          warning("Some x or y were to equal 0 ! Scaling has been adjusted to correct this")
        }
      }
    }
    
    return(loc_proj@coords)
  }
}

##Returns loc_grid of lager domain
###Buffer 
generate_larger_domain <- function(loc,buffer){
  unique_long <-  sort(unique(loc[,1]))
  unique_lat <-  sort(unique(loc[,2]))
  
  h_long <- abs(unique_long[1] - unique_long[2]) 
  h_lat <- abs(unique_lat[1] - unique_lat[2]) 
  
  bbox <- matrix(0,nrow=2,ncol=2)
  
  bbox[1,] <- c(min(unique_long),min(unique_lat)) 
  bbox[2,] <- c(max(unique_long),max(unique_lat)) 
  
  diff_long <- h_long*ceiling((bbox[2,1]-bbox[1,1])/h_long * buffer)
  diff_lat <- h_lat*ceiling((bbox[2,2]-bbox[1,2])/h_lat * buffer)
  
  
  larger_unique_long <- seq(from=bbox[1,1]-diff_long,to=bbox[2,1]+diff_long,by=h_long)
  larger_unique_lat <- seq(from=bbox[1,2]-diff_lat,to=bbox[2,2]+diff_lat,by=h_lat)
  
  return(expand.grid(larger_unique_long,larger_unique_lat))
}

####Here it is supposed the larger and smaller region are squared and ordered as in expand.grid
get_region_indices <- function(large_loc,small_loc){
  bbox <- matrix(0,nrow=2,ncol=2)
  
  bbox[1,] <- c(min(small_loc[,1]),min(small_loc[,2])) 
  bbox[2,] <- c(max(small_loc[,1]),max(small_loc[,2])) 
  
  which( ((large_loc[,1] >= bbox[1,1]) & (large_loc[,1] <= bbox[2,1]) & (large_loc[,2] >= bbox[1,2]) & (large_loc[,2] <= bbox[2,2])) ,arr.ind = TRUE)
}