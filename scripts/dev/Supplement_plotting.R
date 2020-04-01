SEA_noNA <- SEA[complete.cases(SEA), ]             # remove rows of mismatched NA, listwise deletion #
LAND_noNA <- LAND[complete.cases(LAND), ]

str(SEA_submission)
str(LAND_submission)

ggplot(SEA_submission,aes(y = y, x = marine_richness))+
  geom_point(size = .2)+themeo

ggplot(LAND_submission,aes(y = y, x = terrestrial_richness))+
  geom_point(size = .2)+themeo

# %>% slice( sample(1:n()))
a<-ggplot(LAND_submission %>% slice( sample(1:n()))
          ,aes(y = y, x = terrestrial_richness, fill = terrestrial_richness))+
  scale_fill_distiller(palette = "Spectral")+
  geom_point(size = 2.5, shape = 21,show.legend = F)+themeo

b<-ggplot(SEA_submission %>% slice( sample(1:n()))
          ,aes(y = y, x = marine_richness, fill = marine_richness))+
  scale_fill_distiller(palette = "Spectral")+
  geom_point(size = 2.5, shape = 21,show.legend = F)+themeo

c<-ggplot(LAND_submission %>% slice( sample(1:n()))
          ,aes(y = y, x = modeled, fill = modeled))+
  scale_fill_distiller(palette = "Spectral")+
  geom_point(size = 2.5, shape = 21,show.legend = F)+themeo

d<-ggplot(SEA_submission %>% slice( sample(1:n()))
          ,aes(y = y, x = modeled, fill = modeled))+
  scale_fill_distiller(palette = "Spectral")+
  geom_point(size = 2.5, shape = 21,show.legend = F)+themeo
grid.arrange(a,b,c,d)

ggplot(SEA_submission, aes(x = y, y = marine_richness)) + 
  geom_point(size = .2, color = "blue4") +
  geom_line(aes(y=zoo::rollmean(marine_richness, 100, na.pad=TRUE))) +
  themeo

ggplot()+
  geom_point(data = LAND_submission,aes(y = y, x = terrestrial_richness),size = .2, color = "green4")+
  geom_point(data = SEA_submission,aes(y = y, x = marine_richness),size = .2, color = "blue4", alpha = .5)+
  themeo

ggplot()+
  geom_point(data = LAND_submission,aes(y = y, x = terrestrial_richness),size = .2, color = "green4")+
  geom_point(data = SEA_submission,aes(y = y, x = marine_richness),size = .2, color = "blue4", alpha = .5)+
  themeo

ggplot()+
  geom_point(data = LAND_submission,aes(y = y, x = modeled),size = .2, color = "green4")+
  geom_point(data = SEA_submission,aes(y = y, x = modeled),size = .2, color = "blue4", alpha = .5)+
  themeo

ggplot()+
  #geom_point(data = LAND_submission, aes(x = y, y = terrestrial_richness), color = "green4", alpha = .5)+
  #geom_point(data = SEA_submission, aes(x = y, y = marine_richness), color = "blue4", alpha = .5)+
  geom_smooth(data = LAND_submission, aes(x = y, y = terrestrial_richness), color = "green4", alpha = .5)+
  geom_smooth(data = SEA_submission, aes(x = y, y = marine_richness), color = "blue4", alpha = .5)+
  coord_flip()+
  themeo

d_real <- LAND_submission %>%  dplyr::mutate(b=cut(y, breaks=seq(min(LAND_submission$y),max(LAND_submission$y), length.out = 50))) %>%  # 100 cuts is what is Fig 1. explored 50
  dplyr::group_by(b) %>%  dplyr::summarize(median=mean(terrestrial_richness))

m_real <- SEA_submission %>%  dplyr::mutate(b=cut(y, breaks=seq(min(LAND_submission$y),max(LAND_submission$y), length.out = 50))) %>% 
  dplyr::group_by(b) %>%  dplyr::summarize(median=mean(marine_richness))

d_real$b <- as.numeric(d_real$b)
m_real$b <- as.numeric(m_real$b)

m_real$median <- (m_real$median - min(m_real$median))/ (max(m_real$median)-min(m_real$median))
d_real$median <- (d_real$median - min(d_real$median))/ (max(d_real$median)-min(d_real$median))

histo <- ggplot() + 
  geom_bar(data=d_real, aes(x=b, y=median ), stat = 'identity', fill = "green4", width = 1)+
  geom_bar(data=m_real, aes(x=b, y=median ), stat = 'identity', fill = "blue4", width = 1, alpha = .5)+
  theme (axis.text.x=element_text(angle=90, size=5, face='bold'))+
  #theme_void()+
  coord_flip()
histo

mappo <- ggplot()+
  geom_raster(data = SEA_submission, aes(x,y, fill = marine_richness), show.legend = F)+coord_fixed()+
  geom_raster(data = LAND_submission, aes(x,y, fill = terrestrial_richness), show.legend = F)+coord_fixed()+
  # geom_path(data = countries_df,aes(x= long, y = lat, group = group), color = "black", size = .5)+
  scale_fill_distiller(palette = "Spectral", direction = -1)+ theme_void()

grid.arrange(grobs = gList(ggplotGrob(mappo), ggplotGrob(histo)), layout_matrix = cbind(1,1,1,1,1,1,1,1,1,1,1,1,2))

str(d)

mappo

ggplot()+
  geom_raster(data = SEA_submission, aes(x,y, fill = modeled), show.legend = F)+
  geom_raster(data = LAND_submission, aes(x,y, fill = modeled), show.legend = F)+
  geom_path(data = countries_df,aes(x= long, y = lat, group = group), color = "black", size = .25)+
  themeo+
  scale_x_continuous(expand = c(-0.005,-0.005)) +
  scale_y_continuous(expand = c(-0.01,-0.01)) +
  scale_fill_distiller(palette = "Spectral", direction = -1, na.value = "transparent") +
  #theme(axis.text = element_blank(),
  #    axis.title = element_blank())+
  ggtitle("A. Modeled terrestrial and marine richness")+
  coord_fixed()


# modeled richness histo
d <- LAND_submission %>%  dplyr::mutate(b=cut(y, breaks=seq(min(LAND_submission$y),max(LAND_submission$y), length.out = 50))) %>% 
  dplyr::group_by(b) %>%  dplyr::summarize(median=mean(modeled))

m <- SEA_submission %>%  dplyr::mutate(b=cut(y, breaks=seq(min(LAND_submission$y),max(LAND_submission$y), length.out = 50))) %>% 
  dplyr::group_by(b) %>%  dplyr::summarize(median=mean(modeled))

str(d)
str(m)

d$b <- as.numeric(d$b)
m$b <- as.numeric(m$b)

m$median <- (m$median - min(m$median))/ (max(m$median)-min(m$median))
d$median <- (d$median - min(d$median))/ (max(d$median)-min(d$median))

land_histo <- ggplot() + 
  geom_bar(data=d_real, aes(x=b, y=median ), stat = 'identity', fill = "green4", width = 1)+
  geom_bar(data=d, aes(x=b, y=median ), stat = 'identity', fill = "gray", width = 1, alpha = .5)+
  theme (axis.text.x=element_text(angle=90, size=5, face='bold'))+
  theme_void()+
  scale_y_reverse()+
  scale_x_continuous(limits = c(0,50))+
  coord_flip()
land_histo

water_histo <- ggplot() + 
  geom_bar(data=m_real, aes(x=b, y=median), stat = 'identity', fill = "blue4", width = 1)+
  geom_bar(data=m, aes(x=b, y=median), stat = 'identity', fill = "gray", width = 1, alpha = .5)+
  theme (axis.text.x=element_text(angle=90, size=5, face='bold'))+
  theme_void()+
  #scale_y_reverse()+
  scale_x_continuous(limits = c(0,50))+
  coord_flip()
water_histo

grid.arrange(land_histo,water_histo, ncol = 2)






