##  home made functions ##

library(terra)


#######  prepare spatial data given parameters ########

harmonize_raster <- function(x, PU, data_type = c("conti", "categ")){ ## continuous or categorical

  if(data_type == "conti"){
    method_rp <- "bilinear"
    if ( res(x)[1] < res(PU)[1] ){ ## aggregate to coarser resolution if needed before reprojection
      f<- (res(PU)[1])/(res(x)[1])
      x<- terra::aggregate(x, fact = f, fun = "sum", na.rm=T)
      x <- x/(f^2)
    }
  }
  ## x is the raw raster to convert
  ## PU is a raster with the right resolution, projection, extent, origin
  ## Study Area is a vector that has the outline of the spatial area included (e.g. member states)

  ## resolution
  else if (data_type == "categ") {
    method_rp <- "near"
    if ( res(x)[1] < res(PU)[1] ){
      f<- (res(PU)[1])/(res(x)[1])
      x<- terra::aggregate(x, fact = f, fun = "modal", na.rm=T)
    }
  }

  ## project to right grid
  x <- terra::project(x, PU, method = method_rp)

  ### mask out unnecessary countries based on shapefile of study area
  x <- mask(x, PU)

  if (!is.factor(x[cells(x)][1, 1])){
    ### fill in the gaps to have the exact same number of grid cells
    x <- sum(x, PU, na.rm = T)
  }
  return (x)
}



## create a function that is adaptable to different groupings:
summarize_rpz <- function(problem, solutions, feature_info_table, group.by){
  df <- data.frame(solutions = character(),
                   group.by = character(),
                   mean_rpz = numeric())

  for (i in 1:nlyr(solutions)){
    rpz.s_i <- eval_feature_representation_summary(p, solutions[[i]]) ## for each species

    if(group.by ==  "all"){
      df_i <- data.frame(solutions = names(solutions)[i],
                         group.by = group.by,
                         mean_rpz = mean(rpz.s_i$relative_held))

      df <- rbind(df, df_i)
    }
    else{
      rpz.s_i[, group.by] <- feature_info_table[, group.by][match(rpz.s_i$feature, feature_info_table$spp_name)]
      df_i <- aggregate(x = rpz.s_i[c("relative_held")],
                        by = rpz.s_i[c(group.by)],
                        FUN = "mean")

      df_i$solutions <- names(solutions)[i]
      df_i <- df_i[, c(3, 1, 2)]
      colnames(df_i ) <- c("solutions", "group.by", "mean_rpz")
      df <- rbind(df, df_i)

    }
  }

  return(df)
}
