# Vectorized function to split a single geometry matrix 
  split_single_line <- function(line_geom, max_len = 90) {
    
    #get coords
    coords <- st_coordinates(line_geom)[, 1:2]
    
    #calc length
    dists <- c(0, cumsum(sqrt(diff(coords[,1])^2 + diff(coords[,2])^2)))
    total_len <- max(dists)
    
    #if less than max_len, return
    if (total_len <= max_len) return(st_geometry(line_geom))
    
    #break into sections by max_len
    breaks <- seq(0, total_len, by = max_len)
    if (tail(breaks, 1) < total_len) breaks <- c(breaks, total_len)
    
    #find breaks
    new_x <- approx(dists, coords[,1], xout = breaks)$y
    new_y <- approx(dists, coords[,2], xout = breaks)$y
    
    #turn back into a sf object
    map(1:(length(breaks) - 1), ~ {
      st_linestring(matrix(c(new_x[.x:(.x+1)], new_y[.x:(.x+1)]), ncol = 2))
    }) %>% st_sfc(crs = st_crs(line_geom))
  }
