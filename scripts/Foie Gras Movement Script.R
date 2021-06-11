#FoieGras
setwd("C:/R/Humpbackwhalemixing")
library(remotes)
library(tidyverse)
#remotes::install_github("ianjonsen/foieGras")
library(foieGras)
library(sf)
#install.packages("rnaturalearth")
library("rnaturalearth")
#install.packages("rnaturalearthdata")
library("rnaturalearthdata")
library(wesanderson)

####load data
Whales= read.csv(file.choose(), sep = ";",header = T)
View(head(Whales))

###########################convert and combine date time
library(lubridate)
Whales$new =with(Whales, ymd(Date, tz ="GMT") + hm(Time) )
View(head(Whales))

South.africa =Whales[,c(1,8,5,7,6)] #Name which dataset (South Africa or Brazil) and filter required columns 
View(head(South.africa)) 
colnames(South.africa) = c("id","date","lc","lon","lat") #change columns name

#
South.africa= South.africa[ South.africa$id == "RSA-2014-142929",]#run one at the time/subset

#|South.africa$id == "RSA-2019-66181"|South.africa$id == "RSA-2016-131135",]
#South.africa= South.africa[ South.africa$id == "RSA-2016-158097",]
#South.africa= South.africa[ South.africa$id == "RSA-2016-131135",]
#South.africa= South.africa[ South.africa$id == "RSA-2014-142932",]
#South.africa= South.africa[ South.africa$id == "RSA-2014-142931",]


fit <- fit_ssm(South.africa,
               model = "rw",
               vmax = 5,#speed
               ang = c(15, 25),#angle
               distlim = c(2500, 5000),#distance limit 
               time.step = 12#time step
               
)


plot(fit, what = "fitted", type = 1)#temporal plot
plot(fit, what = "fitted", type = 2)#spatial

plot(fit, what = "predicted", type = 1)#temporal plot
plot(fit, what = "predicted", type = 2)#spatial



#fit  movement persistent model from output of fit_ssm              
fmp <- fit %>% 
  grab(what = "fitted", as_sf = FALSE) %>%
  select(id, date, lon, lat) %>%
  fit_mpm(model = "jmpm", verbose = 0) ## turn off parameter trace for tidy output

plot(fmp)#plot movement persistant graph TEMPORAL
fmap(fit, fmp, what = "fitted", 
     crs = "+proj=longlat +datum=WGS84")#SPATIAL projection on map neater


#grab the data, extract for plotting
ten =grab(fit, what = "fitted", as_sf = F)#grab data 
preds <- grab(fmp, what = "fitted", as_sf = FALSE)#ith geometry easier to use

super.data = cbind(ten[,1:6],preds[,1:3])#combine the datasets 
write.csv(super.data,file="RSA142931Output.csv")


#Plot graphs on world map


#scale_fill_gradientn(name="Altitude",colours = cols)+
# guides(fill = guide_colorbar()) +
#scale_alpha(range = c(0.8533891, 10.9017914), guide = "none")

world <- ne_countries(scale = "medium", returnclass = "sf") #define world data
class(world)

ggplot(data = world) +theme_bw()+ geom_sf(color = "grey30", fill = "grey30")+
  coord_sf(crs= "+proj=longlat +datum=WGS84")+
  geom_point(data = super.data[,-c(1:2)], aes(lon, lat, group = id, colour = g),size =2)+
  scale_y_continuous(limits = c(-65,-30))+scale_x_continuous(limits = c(-70,37))+
  scale_colour_gradientn(colours = rev(wes_palette(name = "Zissou1", 
                                                   type = "continuous")), name = expression(gamma[t]), limits = c(0, 1))#















