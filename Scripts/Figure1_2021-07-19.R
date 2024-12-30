######
#figure 1

require(rgeos)
require(sp)
require(rgdal)
require(raster)
require(ggplot)
require(ggpubr)
require(maps)
require(sf)
require(RColorBrewer)
require(cowplot)
require(ggimage)
require(magick)

######
#GIS layers

gis.path <- "C:/Users/tosam/Documents/0_OSU/Dissertation/GIS"
gis.path2 <- "C:/Users/tosam/Documents/0_OSU/Dissertation/HJAndrews/GIS"
nad83z10 <- "+proj=utm +zone=10 +datum=NAD83 +units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"

hjaboundary <- readOGR(gis.path2, "boundary83")
roads <- readOGR(gis.path, "roads_bioarea") #bioarea roads
wb <- readOGR(dsn=paste(gis.path, "nhdwaterbody", sep="/"), layer="nhdwaterbody") #waterbodies for WNF
wb <- spTransform(wb, projection(hjaboundary))

e <- raster(paste(gis.path, "/gi00301_DEM/gi00301.e00", sep=""))
#elevation <- raster(paste(gis.path, "/latlong_bare_earth.tif", sep="")) #lat long projection
#e.crop <- crop(elevation, extent(spTransform(hjaboundary, projection(elevation))))
e.df <- as.data.frame(e, xy=T)
names(e.df)[3] <- "value"

#smgrids <- readOGR(gis.path2, "small_mammal_grids")
smgrids <- readOGR(gis.path2, "small_mammal_grid_corners") #grids as squares


hja <- st_as_sf(hjaboundary)
sm <- st_as_sf(smgrids)

p <- ggplot() +
  #geom_tile(data = e.df, aes(x=x, y=y, fill=value)) +
  geom_raster(data=e.df, aes(x=x, y=y, fill=value)) +
  geom_sf(data=hja, fill=NA, lwd=1, col="black") +
  geom_sf(data=sm, fill="white", lwd=0.5) +
  geom_sf_text(data=sm, aes(label=d______), check_overlap=T) +
  scale_fill_gradient(low="grey10", high="white", name="elevation (m)") +
  coord_sf(xlim=c(558992, 571748), ylim=c(4894224, 4903456)) +
  ylab("") + xlab("") +
  theme_bw(base_size = 20) + theme(legend.position="bottom", legend.key.width = unit(0.6,"in"))
p

save_plot("Figures_Final/Figure1a.png", plot=p)
ggsave(p, filename = "Figures_Final/Figures1a.tiff", width=10, height=8, units="in", dpi=300, compression="lzw")

#########
#oregon map and star for HJA
counties <- readOGR(gis.path, "orcntypoly")
c.sf <- st_as_sf(spTransform(counties, projection(hjaboundary)))

hja.center <- gCentroid(hjaboundary)

star <- image_transparent(image_read("https://upload.wikimedia.org/wikipedia/en/f/f4/Free_Blue_Star.jpg"), 'white')
#star <- image_colorize(star, 100, "#3399FF")
#image_write(star, path = "star.png", format = "png")

p2 <- ggplot() +
  geom_sf(data=c.sf) +
  geom_sf(data=hja) +
#  geom_image(data=data.frame(hja.center), aes(x=x, y=y, image="https://upload.wikimedia.org/wikipedia/en/f/f4/Free_Blue_Star.jpg"), size=0.1) +
  ylab("") + xlab("") +
  theme_bw()

ggsave(p2, filename = "Figures_Final/Figures1b.tiff", width=3, height=2, units="in", dpi=300, compression="lzw")
