# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Attempt to categorized the ASDN points by slope buffer
# 2024-03-17 ---- DMP ----------------------------------------------------------------------------
# Modified 2024-04-19 ANK-------------------------------------------------------

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~ Setup Workspace -------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Clear Memory
remove(list=ls())

#Load libraries of interest
library(tidyverse) #join the cult!
library(raster)
library(sf)
library(whitebox)
library(stars)
library(fasterize)
library(mapview)
library(parallel)
library(ggplot2)
library(ggpubr)
library(fpCompare)

#Define directory
data_dir<-"/Users/ashleigh/Documents/InputData/" 
temp_dir<-"/Users/ashleigh/Documents/Rscratch/" 

#get DEM
DEM <- raster(paste0(data_dir,"DEM/ShambleyFill_dem.tif"))

#get Streams
streams <- st_read(paste0(temp_dir, "streams_weyer.shp"))
st_crs(streams) <- st_crs(DEM)
#streams <- streams %>% 
#  st_transform(pts, crs ='+proj=utm +zone=16 +datum=WGS84 +units=m +no_defs +ellps=WGS84')


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#~~~~~~~ Define Functions ------------------------------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#--------Create XS's -----------------------------------------------------------

#Create function to create cross sections
xs_fun<-function(n, pnts, width=16){
  pnt<-pnts[n,]
  
  #Define flowline segment
  reach<-st_intersection(streams, st_buffer(pnt, dist = 1))
  reach_backup<-reach
  
  #Estimate planar slope
  reach<-st_coordinates(reach)
  reach_slope<-(reach[1,"Y"]-reach[nrow(reach),"Y"])/(reach[1,"X"]-reach[nrow(reach),"X"])
  
  #Estimate inverse slope
  xs_slope <- -1/reach_slope
  
  #Estimate endpoints of XS
  xs_coord <- st_coordinates(pnt)
  xs_coord <-rbind(
    xs_coord, 
    matrix(0, nrow=2, ncol=2)
  )
  xs_coord[2,"X"] <- xs_coord[1,"X"] + width/2*cos(atan(xs_slope))
  xs_coord[2,"Y"] <- xs_coord[1,"Y"] + width/2*sin(atan(xs_slope))
  xs_coord[3,"X"] <- xs_coord[1,"X"] - width/2*cos(atan(xs_slope))
  xs_coord[3,"Y"] <- xs_coord[1,"Y"] - width/2*sin(atan(xs_slope))
  xs_coord<-xs_coord[-1,]
  
  #Create XS
  xs<-xs_coord %>%  
    as_tibble() %>% 
    st_as_sf(coords = c("X","Y")) %>% 
    st_coordinates() %>% 
    st_linestring() %>% 
    st_sfc(.) 
  
  st_crs(xs) <- st_crs(streams) 
  
  xs <- xs %>% 
    st_as_sf() 
  
  #Export XS Shape
  xs
}

#--------Extract elevation data ------------------------------------------------

# Create function to extract elevation data
elevation_fun<-function(n, xs_asdn){
  
  #Call Libraries of interest
  library(tidyverse) #join the cult!
  library(raster)
  library(sf)
  
  #isolate xs
  xs_temp<-xs_asdn[n,]
  
  #convert to raster
  m<-DEM*0+1
  m<-crop(m, st_buffer(xs_temp,250))
  xs_temp<-rasterize(xs_temp, m)
  
  #Add values
  xs_temp<-xs_temp*DEM
  
  #Convert to points
  xs_temp<-rasterToPoints(xs_temp)
  
  #Convert to 2D XS
  xs_temp <- xs_temp %>% 
    #Convert to tibble
    as_tibble() %>% 
    #Estimate distance along transect using distance formula
    mutate(
      dist = ((x[1]-x)^2 + (y[1]-y)^2)^0.5
    ) %>% 
    #Select cols of interest
    dplyr::select(dist, ele = layer)
  
  #Interpolate in the x and y directions
  interp_fun<-approxfun(xs_temp$dist, xs_temp$ele)
  xs_temp<-tibble(
    dist = seq(0, max(xs_temp$dist), 1),
    ele  = interp_fun(seq(0, max(xs_temp$dist), 1)), 
    xs_id = n)
  
  #Export data
  xs_temp
}


#------- Estimate XS area, w/d ratio, and geomorphic unit  ---------------------

# Create function
xs_metrics_fun<-function(n, xs_ele){
  
  #Call Libraries of interest
  library(tidyverse) #join the cult!
  
  #Identify XS of interest
  id<-xs_ele %>% dplyr::select(xs_id) %>% unique()
  id<-id$xs_id[n]
  
  #Isolate data
  xs_ele<- xs_ele %>% filter(xs_id == id)
  
  #Find the top of the left bank
  top_ele_left<-xs_ele %>%
    filter(dist>=2, dist<8) %>%
    filter(ele == max(ele)) %>%
    dplyr::select(ele) %>% pull()
  top_ele_left<-top_ele_left[1]
  
  #Find the top of the right bank
  top_ele_right<-xs_ele %>%
    filter(dist >=8, dist<=14) %>%
    filter(ele == max(ele)) %>%
    dplyr::select(ele) %>% pull()
  top_ele_right<-top_ele_right[1]
  
  #Identify the stream bottom
  invert_ele <- xs_ele %>%
    filter(dist>2, dist<14) %>%
    filter(ele == min(ele)) %>%
    dplyr::select(ele) %>% pull()
  invert_ele <- invert_ele[1]
  
  #Find the distance of the lowest point
  invert_dist<-xs_ele %>%
    filter(dist>2, dist<14) %>%
    filter(ele == invert_ele) %>%
    dplyr::select(dist) %>% pull()
  invert_dist <- invert_dist[1]
  
  #Create function to estimate metrics from XS
  metrics_fun<-function(dz, top_ele_right, top_ele_left){
    #select points which are lower in elevation than the lowest bank top
    #this can be modified if we decide to choose points lower than the highest 
    #bank top, or the closest bank top, instead.
    
    df_left<-xs_ele %>%
      filter(dist>2, dist <8) %>%
      filter(ele < min(top_ele_left, top_ele_right))
    
    df_right<- xs_ele %>%
      filter(dist>=8, dist<14) %>%
      filter(ele < min(top_ele_left, top_ele_right))
    
    #combine the points from the left and right side of the channel
    df <- rbind(df_left, df_right)
    print(df)
    print(min(top_ele_left, top_ele_right))
    
    #Define groups of connected points
    df<-df %>%
      mutate(dx = dist-lag(dist),
             dx = replace_na(dx,0)) %>%
      mutate(group_marker = ifelse(dx==1,0,1)) %>%
      mutate(group = cumsum(group_marker))
    
    #Identify area connected to channel which includes the minimum elevation
    xs_group<-df %>%
      filter(dist==invert_dist) %>%
      dplyr::select(group) %>%
      pull()
    
    #Limit xs to channel which includes the minimum elevation
    df<-df %>% filter(group==xs_group) %>% dplyr::select(-group) %>% dplyr::select(-group_marker) %>%
      dplyr::select(-dx)
    
    #Estimate metrics
    output<-df %>%
      mutate(
        dx_lag = (dist - lag(dist))/2,
        dx_lag = replace_na(dx_lag,0),
        dx_lead = (lead(dist)-dist)/2,
        dx_lead = replace_na(dx_lead,0),
        dx = dx_lag + dx_lead,
        #dy is the distance between the max elevation and every point of elevation
        dy = (dz + invert_ele) - ele,
        print(dy),
        dA = dx*dy) %>%
      summarise(
        dz = dz,
        a_xs = sum(dA),
        d_mean = mean(dy),
        d_max = max(dy),
        w_mean = a_xs/d_mean,
        w_max = max(dist)-min(dist),
        w_d_ratio = w_max/dz)
    
    #Export output
    output
  }
  
  #Apply metrics
  dz<-mean(top_ele_left, top_ele_right)-invert_ele
  metrics<-lapply(X=dz,top_ele_right=top_ele_right, top_ele_left=top_ele_left, FUN = metrics_fun) %>% bind_rows()
  
  #Add Unique ID
  metrics<-metrics %>% mutate(xs_id = id)
  
  #Export metrics
  metrics
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#------------------- Deal with the Points on Streams ---------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Define sampling points along stream
south_pnts <- st_read(paste0(data_dir, "South_XS_points.shp"))
st_write(south_pnts, paste0(temp_dir, "\\south_pnts.shp"), append=F)

central_pnts <- st_read(paste0(data_dir, "Central_XS_points.shp"))
st_write(central_pnts, paste0(temp_dir, "\\central_pnts.shp"), append=F)

north_pnts <- st_read(paste0(data_dir, "North_XS_points.shp"))
st_write(north_pnts, paste0(temp_dir, "\\north_pnts.shp"), append=F)

#Snap sampling points to streams
wbt_snap_pour_points(
  pour_pts = "south_pnts.shp",
  flow_accum = "fac.tif",
  snap_dist = 45,
  output = "snap_south_pnts.shp",
  wd = temp_dir
)

wbt_snap_pour_points(
  pour_pts = "central_pnts.shp",
  flow_accum = "fac.tif",
  snap_dist = 45,
  output = "snap_central_pnts.shp",
  wd = temp_dir
)

wbt_snap_pour_points(
  pour_pts = "north_pnts.shp",
  flow_accum = "fac.tif",
  snap_dist = 45,
  output = "snap_north_pnts.shp",
  wd = temp_dir
)

# Read snapped points into R
south_pnts<-st_read(paste0(temp_dir,"//snap_south_pnts.shp"))
central_pnts<-st_read(paste0(temp_dir,"//snap_central_pnts.shp"))
north_pnts<-st_read(paste0(temp_dir,"//snap_north_pnts.shp"))

#Cut off points that are above the stream network
south_pnts <- south_pnts[-c(0:5),]
central_pnts <- central_pnts
north_pnts <- north_pnts[-c(0:6),]

#Cut off points that are below the stream network
south_pnts <- south_pnts[-c(57:65),]
central_pnts <- central_pnts[-c(89:90),]
north_pnts <- north_pnts[-c(44:47),]

#Plot for funzies
mapview(
  DEM,
  alpha.regions=0.9,
  map.types=c("OpenTopoMap")) +
  mapview(streams, color=c("dark blue"))+
  mapview(north_pnts, col.regions=c("dark orange")) +
  mapview(central_pnts, col.regions=c("green")) +
  mapview(south_pnts, col.regions=c("yellow"))


#Rename rows so that they match xs_ids
rownames(south_pnts) <- 1:nrow(south_pnts)
rownames(central_pnts) <- 1:nrow(central_pnts)
rownames(north_pnts) <- 1:nrow(north_pnts)

#back to shapefiles for arcgis

st_write(south_pnts, paste0(data_dir, "south_pnts.shp"), delete_layer = T)
st_write(central_pnts, paste0(data_dir, "central_pnts.shp"), delete_layer = T)
st_write(north_pnts, paste0(data_dir, "north_pnts.shp"), delete_layer = T)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#-------------------Use the Functions we Made-----------------------------------
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#Apply cross section function
xs_asdn_south <- lapply(seq(1, nrow(south_pnts)), xs_fun, 
                        pnts=south_pnts) %>% bind_rows(.)

xs_asdn_central <- lapply(seq(1, nrow(central_pnts)), xs_fun, 
                          pnts=central_pnts) %>% bind_rows(.)

xs_asdn_north <- lapply(seq(1, nrow(north_pnts)), xs_fun,
                        pnts=north_pnts) %>% bind_rows(.)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Elevation Data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Define number of cores
n_cores<-detectCores()-1

#Create cores
cl<-makeCluster(n_cores)

#Send DEM data to cores
clusterExport(cl, c('DEM'))

#Apply function
xs_ele_south <-parLapply(cl, X=seq(1, nrow(xs_asdn_south)), fun = elevation_fun,
                         xs_asdn=xs_asdn_south) %>% bind_rows()

xs_ele_central <-parLapply(cl, X=seq(1, nrow(xs_asdn_central)), fun = elevation_fun,
                         xs_asdn=xs_asdn_central) %>% bind_rows()

xs_ele_north <-parLapply(cl, X=seq(1, nrow(xs_asdn_north)), fun = elevation_fun,
                         xs_asdn=xs_asdn_north) %>% bind_rows()


#Stop clusters
stopCluster(cl)

#Check that we have cross sections
ggplot(subset(xs_ele_south, xs_id == 36), aes(dist, ele)) +
  ggtitle("Wetland Corridor") + geom_line(color='#1b9e77', size=2) + theme_classic(base_size=20)+  
  scale_x_continuous("Distance (m)", breaks=seq(0,16,4))
  

ggplot(subset(xs_ele_central, xs_id == 16), aes(dist, ele)) +
  ggtitle("Incised Corridor") + geom_line(color='#d95f02', size=2) + theme_classic(base_size=20)+
  scale_x_continuous("Distance (m)", breaks=seq(0,16,4))

ggplot(subset(xs_ele_north, xs_id == 3), aes(dist, ele)) +
  ggtitle("Intact Riparian Corridor") + geom_line(color='#7570b3', size=2) + theme_classic(base_size=20)+
  scale_x_continuous("Distance (m)", breaks=seq(0,16,4))


#~~~~Cross Section Metrics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Define number of cores
n_cores<-detectCores()-1

#Create cores
cl<-makeCluster(n_cores)

#Apply function
xs_metrics_south <-lapply(X=seq(1, nrow(xs_asdn_south)), xs_ele = xs_ele_south, 
                         FUN= xs_metrics_fun) %>% bind_rows()

xs_metrics_central <-lapply(X=seq(1, nrow(xs_asdn_central)), xs_ele = xs_ele_central, 
                          FUN= xs_metrics_fun) %>% bind_rows()

xs_metrics_north <-lapply(X=seq(1, nrow(xs_asdn_north)), xs_ele = xs_ele_north, 
                          FUN= xs_metrics_fun) %>% bind_rows()

#Stop clusters
stopCluster(cl)


#Classify South HGFs based on w_d_ratio and maximum elevation change
xs_wd_south <- xs_metrics_south %>% 
  group_by(xs_id) %>% 
  summarise_all(last) %>% 
  select(xs_id, w_max, w_d_ratio)

xs_final_south <- xs_metrics_south %>% 
  group_by(xs_id) %>% 
  summarise(dy_max = max(d_max),
            dy_mean = mean(d_mean)) %>% 
  left_join(., xs_wd_south) %>% 
  mutate(slope = dy_mean/w_max,
         hgf = ifelse(w_d_ratio <= 15 &
                        dy_max > 0.55 |
                        w_d_ratio < 11.5,
                      1, NA),
         
         hgf = ifelse(w_d_ratio > 25 &
                        dy_max < 0.15 |
                        w_d_ratio > 31.5,
                      3, hgf),
         hgf = ifelse(xs_id== 29, 1, hgf),
         hgf = ifelse(xs_id== 25, 2, hgf), 
         hgf = ifelse(is.na(hgf),
                      2, hgf)) %>% 
        
  
  select(xs_id, hgf, dy_max, w_d_ratio, slope, w_max)

xs_final_south <- xs_final_south[!is.na(xs_final_south$w_d_ratio), ]

xs_final_south %>% 
  group_by(hgf) %>% 
  summarise(slope = mean(slope),
            wd = mean(w_d_ratio),
            dy = mean(dy_max))

#Classify Central HGFs based on w_d_ratio and maximum elevation change
xs_wd_central <- xs_metrics_central %>% 
  group_by(xs_id) %>% 
  summarise_all(last) %>% 
  select(xs_id, w_max, w_d_ratio)

xs_final_central <- xs_metrics_central %>% 
  group_by(xs_id) %>% 
  summarise(dy_max = max(d_max),
            dy_mean = mean(d_mean)) %>% 
  left_join(., xs_wd_central) %>% 
  mutate(slope = dy_mean/w_max,
         hgf = ifelse(w_d_ratio <= 15 &
                        dy_max > 0.55 |
                        w_d_ratio < 11.5,
                      1, NA),
         
         hgf = ifelse(w_d_ratio > 25 &
                        dy_max < 0.15 |
                        w_d_ratio > 31.5,
                      3, hgf),
         
         hgf = ifelse(is.na(hgf),
                      2, hgf)) %>% 
  
  select(xs_id, hgf, dy_max, w_d_ratio, slope, w_max)

xs_final_central <- xs_final_central[!is.na(xs_final_central$w_d_ratio), ]

xs_final_central %>% 
  group_by(hgf) %>% 
  summarise(slope = mean(slope),
            wd = mean(w_d_ratio),
            dy = mean(dy_max))

#Classify North HGFs based on w_d_ratio and maximum elevation change
xs_wd_north <- xs_metrics_north %>% 
  group_by(xs_id) %>% 
  summarise_all(last) %>% 
  select(xs_id, w_max, w_d_ratio)

xs_final_north <- xs_metrics_north %>% 
  group_by(xs_id) %>% 
  summarise(dy_max = max(d_max),
            dy_mean = mean(d_mean)) %>% 
  left_join(., xs_wd_north) %>% 
  mutate(slope = dy_mean/w_max,
         hgf = ifelse(w_d_ratio <= 15 &
                        dy_max > 0.55 |
                        w_d_ratio < 11.5,
                      1, NA),
         
         hgf = ifelse(w_d_ratio > 25 &
                        dy_max < 0.15 |
                        w_d_ratio > 31.5,
                      3, hgf),
         
         hgf = ifelse(is.na(hgf),
                      2, hgf)) %>% 
  
  select(xs_id, hgf, dy_max, w_d_ratio, slope, w_max)

xs_final_north <- xs_final_north[!is.na(xs_final_north$w_d_ratio), ]

xs_final_north %>% 
  group_by(hgf) %>% 
  summarise(slope = mean(slope),
            wd = mean(w_d_ratio),
            dy = mean(dy_max))


# for CENTRAL WATERSHED
#     hgf  slope    wd     dy
#     1    0.0667  9.91   0.519 
#     2    0.0336  20.7   0.242 
#     3    0.0118  80.3   0.083
# Broadly, we are typical from map id 2 until xs id 15 (map id 16) and then wetland from 
# xs_id 16 (map id 15) until xs_id 29 (map id 30)
# and incised starting at xs_38 (map id 39)

# for SOUTH WATERSHED
#     hgf  slope    wd     dy
#     1    0.0658  10.2   0.473 
#     2    0.0306  22.9   0.198 
#     3    0.0112  73.8   0.078
# We are Wet until xs id 11 (map ID 21) and then Typical
# from xs_id 12 (map ID 22) until xs_id 34 (map ID 44)
# Incised from 35 (map ID 45) to 51 (map ID 61)

# for NORTH WATERSHED
#     hgf  slope    wd     dy
#     1    0.0415  12.4   0.323 
#     2    0.0320  23.6   0.259 
#     3    0.0106  62.7   0.068

# We are Wet from xs_id 1 (map ID 15) until xs_id 5 (map ID 19) and then Typical
# until xs_id 23 (map ID 37) and then incised from xs_30 to xs_33 (map ID 44 to 47)
ggplot(subset(xs_ele, xs_id == 28), aes(dist, ele)) + geom_line()

