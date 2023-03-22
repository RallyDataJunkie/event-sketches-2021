library(sf)
library(dplyr)
library(raster)
library(rayshader)
library(ggplot2)
library(ggpubr)
library(pracma)

#geojson_filename = 'data/finland_arctic_2021.geojson'

#route_index = 4 # 3, 4, 5, 6 

#geojson_sf = sf::st_read(geojson_filename)
# Drop the Z dimension
#geojson_sf = st_zm(geojson_sf, drop = TRUE, what = "ZM")

#stage_route_gj = geojsonio::geojson_json(geojson_sf[route_index,]$geometry)


#route_basis_sf = geojson_sf[route_index,]

# Detect the UTM zone as an EPSG code
lonlat2UTMzone = function(lonlat) {
  utm = (floor((lonlat[1] + 180) / 6) %% 60) + 1
  if(lonlat[2] > 0) {
    utm + 32600
  } else{
    utm + 32700
  }
}

get_utm_crs = function(routes){
  
  # Keep track of the original proj4 string
  old_crs = st_crs(geojson_sf[1,])$proj4string
  
  sample_location_x = st_coordinates(st_centroid(routes[1,]))[1]
  sample_location_y = st_coordinates(st_centroid(routes[1,]))[2]
  
  # Generate a new projection in the appropriate UTM zone
  crs_zone = lonlat2UTMzone(c(sample_location_x,
                              sample_location_y))
  
  new_proj4_string = st_crs(crs_zone)$proj4string
  
  new_proj4_string
}

get_utm_projection = function(routes){
  
  new_proj4_string = get_utm_crs(routes)
  # Transform the route to the UTM projection
  utm_routes = st_transform(geojson_sf, crs=new_proj4_string)
  
  utm_routes
  # Or should we returned a named list
  # e.g. including the original projection?
  #list(utm_routes = utm_routes, orig_crs=old_crs)
}

#utm_routes = get_utm_projection(geojson_sf)


#route_basis_utm = utm_routes[route_index,]

# Retrieve elevation raster buffered to 0.5km extent
#stage_coords = as.data.frame(sf::st_coordinates(route_basis_sf))


get_dem = function(stage_coords){
  geoviz::mapzen_dem(stage_coords$Y, stage_coords$X,
                            width_buffer=0.5)
}

#dem <- get_dem(stage_coords)

# Project Raster
get_dem_utm=function(dem,geojson_sf){
  projectRaster(dem, crs = get_utm_crs(geojson_sf))
}
#dem_utm <- get_dem_utm(dem,geojson_sf)

library(rLFT)

stepdist = 10
window = 20
get_route_convexity = function(route_basis_utm){
  bct(route_basis_utm,
      # distance between measurements 
      step = stepdist,
      window = window, ridName = "name")
}
#route_convexity <- get_route_convexity(route_basis_utm)


cornerer = function (df, slight_conv=0.01, closeby=25, large_angle = 20){
  df %>%
    mutate(dirChange = sign(ConvexityIndex) != sign(lag(ConvexityIndex))) %>%
    mutate(straightish =  (abs(ConvexityIndex) < slight_conv)) %>%
    # TO DO PREVIOUS WOULD BE USEFUL IF (stepangle < large_angle)
    mutate(dist =  (lead(MidMeas)-MidMeas)) %>%
    mutate(nearby =  dist < closeby) %>%
    mutate(firstish = !straightish &
             ((nearby & !lag(straightish) & lag(dirChange)) |
                # We don't want the previous node nearby
                (!lag(nearby)) ) & !lag(nearby) )
}

tight_gradient = 0.5

get_route_convexity2 = function(route_convexity, dem_utm){
  route_convexity = cornerer(route_convexity)
  
  routepoints_c = subset(route_convexity,
                         select = c('Midpoint_X', 'Midpoint_Y'))
  
  route_convexity$elevation = raster::extract(dem_utm, routepoints_c)
  route_convexity
}


#route_convexity= get_route_convexity2(route_convexity, dem_utm)

elevation_convexity_plot = function(route_convexity, 
                                    route_basis_sf,
                                    segment_length=1000,
                                    signif_conv_index=0.15){
  min_elevation = max(0,
                      100*round(((min(route_convexity$elevation)/100)-1.6)))
  
  max_elevation = max(route_convexity$elevation)
  max_dist = max(route_convexity$MidMeas)
  
  #Fence post...
  segment_length = ifelse(segment_length<max_dist, segment_length, max_dist)
  
  g = ggplot(route_convexity, aes(x = MidMeas, y=elevation)) +
    geom_line(color='grey') + 
    geom_point(data = tidyr::drop_na(route_convexity[(abs(route_convexity$ConvexityIndex))>signif_conv_index,]), 
               aes(color=ConvexityIndex>0), size=1) 
  
  g = g+ggtitle(route_basis_sf[1]$name)
  g = g+ scale_color_manual(labels = c("left", "right"),
                            values = c("blue", "red"))
  g = g + labs(color = "Significant curvature") 
  g = g+ scale_x_continuous(minor_breaks = seq(segment_length, max_dist, segment_length))
  g
}

library(trajr)
library(dplyr)

create_trj = function(route_basis_utm){
  trj <- TrajFromCoords(as.data.frame(st_coordinates(route_basis_utm)))
  trj  = distinct(trj, x, y, .keep_all = TRUE)
  
  # displacement is a complex number, so we can get the actual distance:
  trj$distance = Mod(trj$displacement) 
  
  # Remove rows with zero distance change
  trj =  trj[!duplicated(trj[,c('x','y')]),]
  
  # Find the accumulated distance at each step
  trj$cum_dist = cumsum(trj$distance)
  
  # Step angle in radians relative to previous
  trj$stepangle = c(0, TrajAngles(trj, compass.direction = NULL) * 180 / pi, NA) 
  
  trj$cumstepangle = cumsum(c(0, TrajAngles(trj, compass.direction = NULL) * 180 / pi, NA))
  
  trj$stepheading = c(TrajAngles(trj, compass.direction = 0)* 180 / pi, NA) 
  
  # Find the gradient of the accumulated angle
  trj$step_gradient = pracma::gradient(trj$cumstepangle, trj$cum_dist)
  
  trj = trj %>%
    mutate(dirChange = lead(sign(step_gradient) != sign(lag(step_gradient))))
  
  slight_gradient = 0.25
  
  trj = trj %>% 
    mutate(straightish =  (abs(step_gradient) < slight_gradient))
  
  # Close distance threshold
  closeby = 20
  
  trj = trj %>%
    mutate(nearby = (distance < closeby) ) %>%
    mutate(significant = abs(step_gradient) > tight_gradient ) %>%
    mutate(flowing = !nearby & !lead(nearby) & 
             !straightish & !significant ) %>%
    mutate(firstish = (significant & !lag(significant)) | (!straightish & 
                                                             ((nearby & !lag(straightish) & lag(dirChange)) |
                                                                (!nearby) ) )) %>%
    mutate(firstish = firstish & !(lag(firstish)))
  
  tight_gradient = 0.5
  
  
  trj = trj %>% 
    mutate(tightens = !firstish & significant & lead(nearby) &
             ((sign(lead(step_gradient))==sign(step_gradient)) & (abs(lead(stepangle)) > abs(step_gradient))))
  
  trj = trj %>% 
    mutate(opens = !firstish & significant & lead(nearby) &
             ((sign(lead(step_gradient))==sign(step_gradient)) & (abs(lead(stepangle)) < abs(step_gradient))))
  
  trj
}

trj_route_plot = function(trj){
  g = ggplot(data=trj,
             aes(x=x, y=y)) + geom_path(color='grey') + coord_sf() +
            xlab('')+ylab('')
  
  g + geom_point( size=0.2, color='blue',
                  data=trj[trj$step_gradient>0.2,]) +
    geom_point(size=0.2,
               data=trj[trj$step_gradient<=-0.2,],
               color='red')
}


trj_segment_plot = function(trj, start, end,
                            x='x', y='y',
                            title='', fix_coords=TRUE,
                            rotate=NULL,
                            gradients=0, gradient_size=2,
                            caption=TRUE) {
  
  # Create the route distance filter limits
  segment_filter = trj$cum_dist >= start &
    trj$cum_dist <= end
  
  # Filter the route
  route_segment = trj[segment_filter,]
  
  if (!plyr::empty(route_segment) & !is.null(rotate)) {
    route_segment = TrajRotate(route_segment,
                               angle = rotate,
                               relative = TRUE)
  }
  
  # Generate the stylised route plot
  g = ggplot(route_segment) +
    geom_path(aes_string(x=x, y=y)) +
    # Add a marker to show the start of the segment
    geom_point(data=head(route_segment,n=1),
               aes_string(x=x, y=y)) +
    theme_void()
  
  # Add a title
  title=as.character(title)
  if (startsWith(title,'auto::')) {
    title = stringr::str_split(title,'::')[[1]][2]
    title_ = paste0(start/1000, '-', end/1000, 'km')
    if (title!='')
      title = paste(title, title_)
    else
      title = title_
  }
  
  if (title!='')
    g = g + ggtitle(title)
  
  # Add a caption
  if (caption) {
    trj_sinuosity = round(TrajSinuosity2(route_segment), digits=2)
    trj_straightness = round(TrajStraightness(route_segment), digits=2)
    g = g + labs(caption = paste("Sinuosity:", trj_sinuosity,
                                 "\nStraight:", trj_straightness))
  }
    
  
  if (fix_coords)
    g = g + coord_fixed()
  
  if (gradients)
    g = g+ geom_point(aes_string(x=x, y=y),
                      size=gradient_size, color='blue',
                      data=route_segment[route_segment$step_gradient>gradients,]) +
    geom_point(aes_string(x=x, y=y), size=gradient_size,
               data=route_segment[route_segment$step_gradient<=-gradients,],
               color='red')
  g
}

trj_segment_multiplot = function(trj, i, title='',
                                 x='x', y='y',
                                 final=FALSE,
                                 step_length = 20,
                                 segment_length=1000,
                                 gradients=0,gradient_size=1,
                                 fix_coords=FALSE, rotate=NULL){
  
  # Preface the start of the stage with a 20m lead
  start_prefix = step_length
  start = segment_length*(i-1)-start_prefix
  if (final) 
    end = Inf
  else
    end = (segment_length*i)
  
  trj_segment_plot(trj, start, end,  x=x, y=y,
                   title=title,
                   fix_coords=fix_coords,
                   rotate=rotate,
                   gradients=gradients, gradient_size=gradient_size)
}

trj_segments_plots =function(trj, segment_length=1000){
  # Create a list to hold each plot as a separate item
  trj_segment_plots = list()
  
  kms = ceil(max(trj$cum_dist)/segment_length)
  # Last kilometer on its own > 350m long
  #if ((max(trj$cum_dist) %% 1000) > 350)
  #  kms = kms+1
  # Iterate through each kilometer
  for (i in 1:kms){
    # Add each plot to the plot list
    trj_segment_plots[[length(trj_segment_plots) + 1]] <-
      trj_segment_multiplot(trj, i,
                            title=i,
                            #final=(i==kms),
                            fix_coords=TRUE, rotate=0,
                            gradients=0.2, gradient_size=1)
  }
  trj_segment_plots
}


#http://michaeljw.com/blog/post/subchunkify/
#generates a new chunk with the desired width and height
subchunkify <- function(g, fig_height=7, fig_width=5, prefix="") {
  g_deparsed <- paste0(deparse(
    function() {g}
  ), collapse = '')
  
  sub_chunk <- paste0("
  `","``{r sub_chunk_",prefix, floor(runif(1) * 1000000), ", fig.height=", fig_height, ", fig.width=", fig_width, ", echo=FALSE}",
                      "\n(", 
                      g_deparsed
                      , ")()",
                      "\n`","``
  ")
  
  cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
}

trj_segments_plot = function(trj_segment_plots, route_basis_sf, ncol=5){
  
  gg = ggarrange(plotlist=trj_segment_plots,
                 ncol=ncol, nrow=ceiling(length(trj_segment_plots)/ncol))
  
  gg = annotate_figure(gg,
                  top = text_grob(route_basis_sf[1]$name, color = "black",
                                  face = "bold", size = 14))

  gg
}


trj_sf = function(trj){
  trj %>% sf::st_as_sf(coords = c("x","y")) %>% 
    sf::st_set_crs(st_crs(utm_routes[route_index,]))
}


library(amt)

to_amt_track = function(route){
  make_track(st_coordinates(route$geometry), X, Y)
}

get_amt_stats = function(trj_km) {
  amt_tracks = apply(trj_km, 1, to_amt_track)
  
  amt_stats = data.frame(
    km = 1:length(amt_tracks),
    amt_sin = unlist(lapply(amt_tracks, amt::sinuosity)),
    amt_str = unlist(lapply(amt_tracks, amt::straightness)),
    amt_cumd = unlist(lapply(amt_tracks, amt::cum_dist)),
    amt_totd = unlist(lapply(amt_tracks, amt::tot_dist)),
    amt_int = unlist(lapply(amt_tracks, amt::intensity_use))
  )
  
  amt_stats
}


curvature = function(x,y){
  #729181.8, 729186.1, 729190.4
  #4957667 , 4957676, 4957685
  tryCatch({
    # circlefit gives an error if we pass a straight line
    # Also hide the print statement in circlefit
    # circlefit() returns the x and y coords of the circle center
    # as well as the radius of curvature
    # We could then also calculate the angle and arc length
    pracma::circlefit(x,y)[3]
  },
  error = function(err) { 
    # For a straight, return the first co-ord and Inf diameter
    # Alternatively, pass zero diameter?
    c(x[1], y[1], Inf)[3]})
}

curvature2 = function(x1, x2, x3, y1, y2, y3){
  curvature(c(x1, x2, x3), c(y1, y2, y3))
}

# The base::Vectorize function provides a lazy way of 
# vectorising a non-vectorised function
curvatures = Vectorize(curvature2)


library(formattable)

start_location_sf = function(route) {
  as.data.frame(st_coordinates(route$geometry))[1,]
}

# Simpler chart
seg_gg = function(route){
  ggplot(st_linestring(cbind(x=route$x,
                             y=route$y))) +
    # Add a marker to show start of section
    geom_point(data=route[1,], aes(x=x, y=y)) +
    geom_sf() +
    theme_void()
}
